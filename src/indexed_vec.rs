//! Type-safe indexing for cells, vertices, boxes and edges.  Each type of object gets its own
//! index type, which is an opaque new-type over integers.  Additionally, each type of object gets
//! its own replacement for [`Vec`], which can only be indexed by its corresponding index type.

use std::{
    fmt::{Debug, Formatter},
    iter::FromIterator,
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

#[allow(dead_code)] // Some functions exist solely to make the API more complete and may be unused
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

    /// Creates a `TypedVec` made of repeating the same element to a given length.  The untyped
    /// equivalent is `vec![elem; len]`.
    pub fn repeat(elem: T, len: usize) -> Self
    where
        T: Clone,
    {
        Self {
            inner: vec![elem; len],
            _phantom_data: PhantomData,
        }
    }

    /// Returns `true` if this `TypedVec` has no elements
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Returns the number of elements in this `TypedVec`
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

    /// Get the element at a given index, or `None` if the index doesn't exist.
    pub fn get(&self, idx: Idx) -> Option<&T>
    where
        Idx: IdxType,
    {
        self.inner.get(idx.to_idx())
    }

    /// Gets a mutable reference to the element at a given index, or `None` if the index doesn't
    /// exist.
    pub fn get_mut(&mut self, idx: Idx) -> Option<&mut T>
    where
        Idx: IdxType,
    {
        self.inner.get_mut(idx.to_idx())
    }

    /* ITER FUNCTIONS */

    /// Returns an [`Iterator`] which yields references to the elements within this `TypedVec`
    pub fn iter(&self) -> std::slice::Iter<T> {
        self.inner.iter()
    }

    /// Returns an [`Iterator`] which yields mutable references to the elements within this
    /// `TypedVec`
    pub fn iter_mut(&mut self) -> std::slice::IterMut<T> {
        self.inner.iter_mut()
    }

    /// Consumes this `TypedVec` and turns it into an [`Iterator`] which yields all the contained
    /// values
    pub fn into_iter(self) -> std::vec::IntoIter<T> {
        self.inner.into_iter()
    }

    /// Returns an [`Iterator`] which yields references to the elements within this `TypedVec`,
    /// along with their indices.  The untyped equivalent is `.iter().enumerate()`.
    pub fn indexed_iter(&self) -> impl Iterator<Item = (Idx, &T)>
    where
        Idx: IdxType,
    {
        self.inner
            .iter()
            .enumerate()
            .map(|(i, v)| (Idx::from_idx(i), v))
    }

    /// Returns an [`Iterator`] which yields mutable references to the elements within this
    /// `TypedVec`, along with their indices.  The untyped equivalent is `.iter().enumerate()`.
    pub fn indexed_iter_mut(&mut self) -> impl Iterator<Item = (Idx, &mut T)>
    where
        Idx: IdxType,
    {
        self.inner
            .iter_mut()
            .enumerate()
            .map(|(i, v)| (Idx::from_idx(i), v))
    }

    /// Returns a new `TypedVec` where a given function has been applied to each element in `self`.
    /// The untyped equivalent is `.iter().map(f).collect::<Vec<_>>()`.
    pub fn map<U>(&self, f: impl Fn(&T) -> U) -> TypedVec<Idx, U> {
        TypedVec {
            inner: self.inner.iter().map(f).collect_vec(),
            _phantom_data: PhantomData,
        }
    }

    /// Gets the index of the first element `x` for which `f(x)` returns `true`.
    pub fn idx_of_first(&self, f: impl Fn(&T) -> bool) -> Option<Idx>
    where
        Idx: IdxType,
    {
        let pos = self.iter().position(f);
        pos.map(Idx::from_idx)
    }
}

impl<IdxT: IdxType, T> From<Vec<T>> for TypedVec<IdxT, T> {
    fn from(contents: Vec<T>) -> Self {
        Self {
            inner: contents,
            _phantom_data: PhantomData,
        }
    }
}

impl<IdxT: IdxType, T> FromIterator<T> for TypedVec<IdxT, T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self {
            inner: iter.into_iter().collect_vec(),
            _phantom_data: PhantomData,
        }
    }
}

impl<IdxT: IdxType, T> Index<IdxT> for TypedVec<IdxT, T> {
    type Output = T;

    fn index(&self, index: IdxT) -> &Self::Output {
        &self.inner[index.to_idx()]
    }
}

impl<IdxT: IdxType, T> IndexMut<IdxT> for TypedVec<IdxT, T> {
    fn index_mut(&mut self, index: IdxT) -> &mut Self::Output {
        &mut self.inner[index.to_idx()]
    }
}

///////////////////////////
// MACRO/TRAIT MACHINERY //
///////////////////////////

macro_rules! idx_impl {
    ($idx_name: ident, $vec_name: ident) => {
        // Define index type
        #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $idx_name {
            idx: usize,
        }

        impl IdxType for $idx_name {
            fn from_idx(idx: usize) -> Self {
                Self { idx }
            }

            fn to_idx(self) -> usize {
                self.idx
            }
        }

        impl Debug for $idx_name {
            fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
                write!(f, "{}({})", stringify!($idx_name), self.idx)
            }
        }

        // Define vec type alias
        pub type $vec_name<T> = TypedVec<$idx_name, T>;
    };
}

idx_impl!(BoxIdx, BoxVec); // Collection of boxes of cells
idx_impl!(CellIdx, CellVec); // Collection of cells
idx_impl!(EdgeIdx, EdgeVec); // Collection of edges between either boxes or cells
idx_impl!(VertIdx, VertVec); // Collection of vertices of either boxes or cells (or both)
idx_impl!(SymmIdx, SymmVec); // Collection of equivalence classes of anything
idx_impl!(LinkIdx, LinkVec); // Collection of links between non-adjacent edges

/// A common trait implemented by all custom index types
pub trait IdxType {
    fn from_idx(idx: usize) -> Self;

    fn to_idx(self) -> usize;
}
