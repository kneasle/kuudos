//! Type-safe indexing for cells, vertices, boxes and edges.  Each type of object gets its own
//! index type, which is an opaque new-type over integers.  Additionally, each type of object gets
//! its own replacement for [`Vec`], which can only be indexed by its corresponding index type.

use std::{
    marker::PhantomData,
    ops::{Index, IndexMut},
};

use itertools::Itertools;

/// A new-type over [`Vec`] which will only accept indices of an opaque index type
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct TypedVec<Idx, T> {
    inner: Vec<T>,
    _phantom_data: PhantomData<Idx>,
}

impl<Idx, T> TypedVec<Idx, T> {
    /// Creates a new, empty type-safe collection
    pub fn new() -> Self {
        Self {
            inner: Vec::new(),
            _phantom_data: PhantomData,
        }
    }

    /// Creates a new, empty type-safe collection which can take `cap` items without reallocating.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            inner: Vec::with_capacity(cap),
            _phantom_data: PhantomData,
        }
    }

    pub fn repeat(elem: T, len: usize) -> Self
    where
        T: Clone,
    {
        Self {
            inner: vec![elem; len],
            _phantom_data: PhantomData,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Gets the index of the next element to be [`push`](Self::push)ed to this collection.
    pub fn next_idx(&self) -> Idx
    where
        Idx: IdxType,
    {
        Idx::from_idx(self.inner.len())
    }

    /// Adds a new element to this collection, returning its opaque index
    pub fn push(&mut self, t: T) -> Idx
    where
        Idx: IdxType,
    {
        let idx = self.next_idx();
        self.inner.push(t);
        idx
    }

    pub fn get(&self, idx: Idx) -> Option<&T>
    where
        Idx: IdxType,
    {
        self.inner.get(idx.as_idx())
    }

    pub fn get_mut(&mut self, idx: Idx) -> Option<&mut T>
    where
        Idx: IdxType,
    {
        self.inner.get_mut(idx.as_idx())
    }

    /* ITER FUNCTIONS */

    pub fn iter(&self) -> std::slice::Iter<T> {
        self.inner.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<T> {
        self.inner.iter_mut()
    }

    pub fn indexed_iter(&self) -> impl Iterator<Item = (Idx, &T)>
    where
        Idx: IdxType,
    {
        self.inner
            .iter()
            .enumerate()
            .map(|(i, v)| (Idx::from_idx(i), v))
    }

    pub fn map<U>(&self, f: impl Fn(&T) -> U) -> TypedVec<Idx, U> {
        TypedVec {
            inner: self.inner.iter().map(f).collect_vec(),
            _phantom_data: PhantomData,
        }
    }

    /// Gets the index of the first element for which `f` returns `true`
    pub fn idx_of_first(&self, f: impl Fn(&T) -> bool) -> Option<Idx>
    where
        Idx: IdxType,
    {
        let pos = self.iter().position(f);
        pos.map(Idx::from_idx)
    }
}

impl<IdxT: IdxType, T> Index<IdxT> for TypedVec<IdxT, T> {
    type Output = T;

    fn index(&self, index: IdxT) -> &Self::Output {
        self.get(index).unwrap()
    }
}

impl<IdxT: IdxType, T> IndexMut<IdxT> for TypedVec<IdxT, T> {
    fn index_mut(&mut self, index: IdxT) -> &mut Self::Output {
        self.get_mut(index).unwrap()
    }
}

///////////////////////////
// MACRO/TRAIT MACHINERY //
///////////////////////////

macro_rules! idx_impl {
    ($idx_name: ident, $vec_name: ident) => {
        /// An index type used for referring to vertices
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $idx_name {
            idx: usize,
        }

        impl IdxType for $idx_name {
            fn from_idx(idx: usize) -> Self {
                Self { idx }
            }

            fn as_idx(self) -> usize {
                self.idx
            }
        }

        pub type $vec_name<T> = TypedVec<$idx_name, T>;
    };
}

idx_impl!(BoxIdx, BoxVec);
idx_impl!(CellIdx, CellVec);
idx_impl!(EdgeIdx, EdgeVec);
idx_impl!(VertIdx, VertVec);

/// A common trait implemented by all custom index types
pub trait IdxType {
    fn from_idx(idx: usize) -> Self;

    fn as_idx(self) -> usize;
}
