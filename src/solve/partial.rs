use std::{
    fmt::{Display, Formatter},
    iter::FusedIterator,
};

use itertools::Itertools;

use crate::{
    indexed_vec::{CellIdx, CellVec, IdxType},
    shape::Shape,
};

use super::Error;

/// The state of a `Partial`ly solved sudoku grid
#[derive(Debug, Clone)]
pub(super) struct Partial {
    /// Which cells have 'penned' values
    penned_cells: CellVec<Option<usize>>,
    /// How many cells are 'unpenned'.  A `Partial` is solved if this is zero.
    ///
    /// Invariant: `num_unpenned_cells = penned_cells.iter().filter(|c| c.is_none()).count()`
    num_unpenned_cells: usize,

    /// A bitmask for each cell, where a `1` at position `i` means that digit `i` can be placed in
    /// that cell
    pencil_masks: CellVec<usize>,
}

impl Partial {
    /// Create a `Partial` representing a completely empty grid
    pub fn empty(shape: &Shape) -> Partial {
        // `all_options` has `shape.num_symbols` 1s in its lowest bits, and 0s elsewhere
        let all_options = (1 << shape.num_symbols) - 1;
        Self {
            penned_cells: CellVec::repeat(None, shape.num_cells()),
            num_unpenned_cells: shape.num_cells(),
            pencil_masks: CellVec::repeat(all_options, shape.num_cells()),
        }
    }

    /// Create a `Partial` with a given set of clues penned in.
    pub fn from_clues(
        table: &Table,
        shape: &Shape,
        clues: &[Option<usize>],
    ) -> Result<Self, Error> {
        // Check that there are as many clues as cells in the `shape`
        if shape.num_cells() != clues.len() {
            return Err(Error::ClueLenMismatch {
                clue_len: clues.len(),
                num_cells: shape.num_cells(),
            });
        }

        let mut partial = Self::empty(shape);
        for (i, clue) in clues.iter().enumerate() {
            if let Some(digit) = clue {
                partial.pen(table, CellIdx::from_idx(i), *digit);
            }
        }
        Ok(partial)
    }

    /// Solve in-place as many naked singles as possible - i.e. if a cell has only one possible
    /// option, then fill it in.  This returns the number of cells solved.
    pub fn solve_naked_singles(&mut self, table: &Table) -> usize {
        let mut num_digits_penned = 0;
        loop {
            let mut has_naked_singles = false;
            for (cell_idx, _cell) in table.affected_cells.indexed_iter() {
                // If this cell is a naked single then we pen its value
                if let Some(only_digit) = self.naked_single_digit(cell_idx) {
                    has_naked_singles = true;
                    self.pen(table, cell_idx, only_digit);
                    num_digits_penned += 1;
                    // Whenever we pen a cell, check that the invariants are still upheld (but only
                    // in debug mode)
                    self.debug_assert_invariants();
                }
            }
            // If we looked at all cells without solving any naked singles, then no more naked
            // singles are possible
            if !has_naked_singles {
                return num_digits_penned;
            }
        }
    }

    /// Pen a value in a given cell
    pub fn pen(&mut self, table: &Table, cell_idx: CellIdx, digit: usize) {
        // Record the penned value in the current cell
        if self.penned_cells[cell_idx].is_none() {
            self.num_unpenned_cells -= 1;
        }
        self.penned_cells[cell_idx] = Some(digit);
        // Set all bits of the pencil mask, to prevent this cell from being pencilled again
        self.pencil_masks[cell_idx] = usize::MAX;
        // 'Unpencil' this value from any cell which shares a house with this one
        for &affected_cell_idx in &table.affected_cells[cell_idx] {
            let pencil_mask = &mut self.pencil_masks[affected_cell_idx];
            let bit = 1 << digit;
            // Set the bit corresponding to this value to 0
            *pencil_mask &= !bit;
        }
    }

    /// Is the given cell a 'naked single' (_\*sniggering intensifies\*_) - i.e. is there only one
    /// cell that can go in it?
    pub fn naked_single_digit(&self, cell_idx: CellIdx) -> Option<usize> {
        let pencil_mask = self.pencil_masks[cell_idx];
        if pencil_mask.count_ones() == 1 {
            Some(63 - pencil_mask.leading_zeros() as usize)
        } else {
            None
        }
    }

    /// Is this grid solved?
    pub fn is_solved(&self) -> bool {
        self.num_unpenned_cells == 0
    }

    /// Returns the digits in this `Partial`, returning `None` if the grid isn't solved.
    pub fn into_solved_digits(self) -> Option<Vec<usize>> {
        self.penned_cells.into_iter().collect()
    }

    /// Gets the [`CellIdx`] of a cell with a minimal number of pencilled digits.
    pub fn least_pencilled_cell(&self) -> CellIdx {
        self.pencil_masks
            .indexed_iter()
            .min_by_key(|(_cell_idx, pencil_mask)| pencil_mask.count_ones())
            .map(|(cell_idx, _)| cell_idx)
            .expect("Grid has no cells")
    }

    /// Gets the pencil mask for a given cell - i.e. a number where the locations of `1`s show
    /// which digits are possible in this cell.
    pub fn pencilled_digit_iter(&self, cell_idx: CellIdx) -> Option<PencilledDigitIter> {
        let pencil_mask = *self.pencil_masks.get(cell_idx)?;
        Some(PencilledDigitIter { pencil_mask })
    }

    /// Asserts that `self` upholds the required invariants, and panics otherwise (does nothing in
    /// release mode)
    #[cfg(debug_assertions)]
    pub fn debug_assert_invariants(&self) {
        assert_eq!(self.pencil_masks.len(), self.penned_cells.len());

        assert_eq!(
            self.num_unpenned_cells,
            self.penned_cells.iter().filter(|c| c.is_none()).count()
        );
    }

    /// Asserts that `self` upholds the required invariants, and panics otherwise (does nothing in
    /// release mode)
    #[cfg(not(debug_assertions))]
    #[inline(always)]
    pub fn debug_assert_invariants(&self) {}
}

impl Display for Partial {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (pen, pencil_mask) in self.penned_cells.iter().zip_eq(self.pencil_masks.iter()) {
            match pen {
                Some(p) => writeln!(f, "{:>2}: {:0>16b}", p, pencil_mask)?,
                None => writeln!(f, "  : {:0>16b}", pencil_mask)?,
            }
        }
        write!(f, "{} cells unpenned", self.num_unpenned_cells)
    }
}

/// Lookup tables used by [`Partial`]s
#[derive(Debug, Clone)]
pub(super) struct Table {
    /// For a cell with index `i`, `affected_cells[i]` contains the set of cells which share a
    /// group with cell `i` (**excluding** `i` itself)
    affected_cells: CellVec<Vec<CellIdx>>,
}

impl Table {
    pub fn from_shape(shape: &Shape) -> Self {
        let mut affected_cells: CellVec<Vec<CellIdx>> =
            CellVec::repeat(Vec::new(), shape.num_cells());
        for group in &shape.groups {
            for &cell_idx in &group.cells {
                affected_cells[cell_idx].extend(group.cells.iter().filter(|idx| **idx != cell_idx));
            }
        }
        // Deduplicate the lists of affected cells (because pairs of cells often have several
        // groups in common)
        for cells in affected_cells.iter_mut() {
            cells.sort_unstable();
            cells.dedup();
        }
        Self { affected_cells }
    }
}

/// An [`Iterator`] which yields the pencilled digits in a given cell (in increasing order)
#[derive(Debug, Clone, Copy)]
pub struct PencilledDigitIter {
    /// Mask where every the index `1` is a pencilled digit which hasn't yet been returned
    pencil_mask: usize,
}

impl Iterator for PencilledDigitIter {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pencil_mask == 0 {
            return None;
        }
        let next_digit = self.pencil_mask.trailing_zeros() as usize;
        self.pencil_mask &= !(1 << next_digit); // Set this digit's bit to `0`
        Some(next_digit)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let digits_left = self.len();
        (digits_left, Some(digits_left))
    }
}

impl ExactSizeIterator for PencilledDigitIter {
    fn len(&self) -> usize {
        self.pencil_mask.count_ones() as usize
    }
}

impl FusedIterator for PencilledDigitIter {}
