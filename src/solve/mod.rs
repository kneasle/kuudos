use std::fmt::{Display, Formatter};

use itertools::Itertools;

use crate::Shape;

/// Solves a sudoku, returning a possible error if there are no or multiple solutions
pub fn solve(
    shape: &Shape,
    clues: &[Option<usize>],
    check_multiple_solns: bool,
) -> Result<Vec<usize>, Error> {
    if shape.num_cells() != clues.len() {
        return Err(Error::ClueLenMismatch {
            clue_len: clues.len(),
            num_cells: shape.num_cells(),
        });
    }

    let data = ShapeData::from(shape);

    // Create a partial solution where only the given clues are penned
    let mut partial = Partial::empty(shape);
    clues
        .iter()
        .enumerate()
        .filter_map(|(i, maybe_v)| maybe_v.map(|v| (i, v)))
        .for_each(|(cell_idx, value)| partial.pen(&data, cell_idx, value));

    // Run recursive backtracking on this grid, and extract the solution or propagate the error
    recursive_solve(&data, partial, check_multiple_solns).map(Partial::into_solved_digits)
}

/// A very naive backtracking solver, which looks for 'naked singles'
fn recursive_solve(
    data: &ShapeData,
    mut partial: Partial,
    check_multiple_solns: bool,
) -> Result<Partial, Error> {
    // Repeatedly pen in 'naked singles' - i.e. a cell with only one possible digit.  These are
    // cheap to compute and don't require backtracking, so we perform them in-place on the current
    // Partial.
    loop {
        let mut has_naked_singles = false;
        for cell_idx in 0..data.num_cells() {
            // If this cell is a naked single then we pen its value
            if let Some(only_digit) = partial.naked_single_digit(cell_idx) {
                has_naked_singles = true;
                partial.pen(data, cell_idx, only_digit);
                // Whenever we pen a cell, check that the invariants are still upheld (but only in
                // debug mode)
                partial.debug_assert_invariants();
            }
        }
        // Exit the loop when no more naked singles are possible
        if !has_naked_singles {
            break;
        }
    }

    // If the grid has not been solved by naked singles, then it is either solved or a branch is
    // required

    // If the value of `partial` passed to this function is soluble using only naked singles then
    // it permits precisely one solution, which `partial` now contains.  So we return that single
    // solution
    if partial.is_solved() {
        return Ok(partial);
    }

    // If `partial` isn't solved but no more naked singles are found, then we pick a non-penned
    // cell and speculatively try each value to see which solutions are found.  We always want to
    // remove as much of the possible search space as possible, so we pick (one of) the cells with
    // the fewest options
    let min_idx = partial
        .num_pencil_marks
        .iter()
        .enumerate()
        .min_by_key(|(_cell_idx, num_opts)| **num_opts)
        .map(|(cell_idx, _)| cell_idx)
        .expect("Grid has no cells");

    let mut solution: Option<Partial> = None;
    let mut unexplored_pencil_mask = partial.pencil_masks[min_idx];
    while unexplored_pencil_mask != 0 {
        // Repeatedly pick the highest pencil mark in the cell, and remove its flag bit
        let assumed_digit = 63 - unexplored_pencil_mask.leading_zeros() as usize;
        unexplored_pencil_mask &= !(1 << assumed_digit);
        // Create a new `Partial` which assumes that this digit is taken
        let mut new_partial = partial.clone();
        new_partial.pen(data, min_idx, assumed_digit);
        match recursive_solve(data, new_partial, check_multiple_solns) {
            Ok(new_solution) => {
                if check_multiple_solns {
                    // If we are required to check for multiple solutions and this is the first
                    // solution found, we continue the search
                    match solution {
                        // If this is the 2nd solution, then this branch is unsolvable
                        Some(_) => return Err(Error::MultipleSolutions),
                        None => solution = Some(new_solution),
                    }
                } else {
                    return Ok(new_solution);
                }
            }
            // If the solver produced too many solutions, then this whole branch is unsolvable
            // (because the other branches can only generate more solutions)
            Err(Error::MultipleSolutions) => return Err(Error::MultipleSolutions),
            // If no solutions where found, then we keep looking for a solution
            Err(Error::NoSolutions) => {}

            // This error can only be generated by the top-level `solve()` function
            Err(Error::ClueLenMismatch { .. }) => unreachable!(),
        }
    }

    // If multiple solutions exist down this branch, we would have taken an early return.
    // Therefore, this branch contains either 0 or 1 solutions (encoded by `solution` being `None`
    // or `Some(_)` respectively).
    solution.ok_or(Error::NoSolutions)
}

/// Useful lookup tables generated from a given [`Shape`]
#[derive(Debug, Clone)]
struct ShapeData {
    /// For a cell with index `i`, `affected_cells[i]` contains the set of cells which share a
    /// house with cell `i` (**excluding** `i` itself)
    affected_cells: Vec<Vec<usize>>,
}

impl ShapeData {
    fn num_cells(&self) -> usize {
        self.affected_cells.len()
    }
}

impl From<&Shape> for ShapeData {
    fn from(shape: &Shape) -> Self {
        let mut affected_cells: Vec<Vec<usize>> = vec![Vec::new(); shape.num_cells()];
        for cell_idxs_in_group in &shape.groups {
            for &cell_idx in cell_idxs_in_group {
                affected_cells[cell_idx]
                    .extend(cell_idxs_in_group.iter().filter(|idx| **idx != cell_idx));
            }
        }
        // Deduplicate the lists of affected cells
        for cells in &mut affected_cells {
            cells.sort_unstable();
            cells.dedup();
        }
        Self { affected_cells }
    }
}

/// The state of a `Partial`ly solved sudoku grid
#[derive(Debug, Clone)]
struct Partial {
    /// Which cells have 'penned' values
    penned_cells: Vec<Option<usize>>,
    /// How many cells are 'unpenned'.
    ///
    /// Invariant: `num_unpenned_cells = penned_cells.iter().filter(|c| c.is_none()).count()`
    num_unpenned_cells: usize,

    /// A bitmask for each cell, where a `1` at position `i` means that digit `i` can be placed in
    /// that cell
    pencil_masks: Vec<usize>,
    /// How many options are left per cell.
    ///
    /// Invariant: for all i: `num_options_left[i] = options_per_cell[i].count_ones()`
    num_pencil_marks: Vec<usize>,
}

impl Partial {
    /// A `Partial` representing a completely empty grid
    fn empty(shape: &Shape) -> Partial {
        // `all_options` has `shape.num_symbols` 1s in its lowest bits, and 0s elsewhere
        let all_options = (1 << shape.num_symbols) - 1;
        Self {
            penned_cells: vec![None; shape.num_cells()],
            num_unpenned_cells: shape.num_cells(),
            pencil_masks: vec![all_options; shape.num_cells()],
            num_pencil_marks: vec![shape.num_symbols; shape.num_cells()],
        }
    }

    /// Pen a value in a given cell
    fn pen(&mut self, data: &ShapeData, cell_idx: usize, value: usize) {
        // Record the penned value in the current cell
        if self.penned_cells[cell_idx].is_none() {
            self.num_unpenned_cells -= 1;
        }
        self.penned_cells[cell_idx] = Some(value);
        // Set the number of pencil marks in this cell to be as large as possible (so this cell
        // will never get chosen for either a hidden single or as the best branch).
        self.num_pencil_marks[cell_idx] = usize::MAX;
        // 'Unpencil' this value from any cell which shares a house with this one
        for &affected_cell_idx in &data.affected_cells[cell_idx] {
            let pencil_mask = &mut self.pencil_masks[affected_cell_idx];
            let bit = 1 << value;
            // If the cell is currently pencilled, then update the corresponding count
            if *pencil_mask & bit != 0 {
                self.num_pencil_marks[affected_cell_idx] -= 1;
            }
            // Set the bit corresponding to this value to 0
            *pencil_mask &= !bit;
        }
    }

    /// Is the given cell a 'naked single' (_\*sniggering intensifies\*_) - i.e. is there only one
    /// cell that can go in it?
    fn naked_single_digit(&self, cell_idx: usize) -> Option<usize> {
        let is_naked_single = self.num_pencil_marks[cell_idx] == 1;
        if is_naked_single {
            Some(63 - self.pencil_masks[cell_idx].leading_zeros() as usize)
        } else {
            None
        }
    }

    /// Is this grid solved?
    fn is_solved(&self) -> bool {
        self.num_unpenned_cells == 0
    }

    /// Returns the digits in this `Partial`, panicking if the grid isn't solved.
    ///
    /// # Panics
    ///
    /// Panics if `self.is_solved()` is `false`
    fn into_solved_digits(self) -> Vec<usize> {
        self.penned_cells
            .into_iter()
            .map(|v| v.expect("`recursive_solve` returned an unsolved grid"))
            .collect_vec()
    }

    /// Asserts that `self` upholds the required invariants, and panics otherwise (does nothing in
    /// release mode)
    #[cfg(debug_assertions)]
    fn debug_assert_invariants(&self) {
        assert_eq!(self.pencil_masks.len(), self.penned_cells.len());
        assert_eq!(self.pencil_masks.len(), self.num_pencil_marks.len());

        assert_eq!(
            self.num_unpenned_cells,
            self.penned_cells.iter().filter(|c| c.is_none()).count()
        );
        for (cell_idx, (pencil_mask, count)) in self
            .pencil_masks
            .iter()
            .zip_eq(&self.num_pencil_marks)
            .enumerate()
        {
            // This invariant only applies to unpenned cells
            if self.penned_cells[cell_idx].is_some() {
                continue;
            }
            assert_eq!(*count, pencil_mask.count_ones() as usize);
        }
    }

    /// Asserts that `self` upholds the required invariants, and panics otherwise (does nothing in
    /// release mode)
    #[cfg(not(debug_assertions))]
    #[inline(always)]
    fn debug_assert_invariants(&self) {}
}

impl Display for Partial {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (pen, (pencil_mask, opts_left)) in self
            .penned_cells
            .iter()
            .zip_eq(self.pencil_masks.iter().zip_eq(&self.num_pencil_marks))
        {
            match pen {
                Some(p) => writeln!(f, "{:>2}: {:0>16b} {}", p, pencil_mask, opts_left)?,
                None => writeln!(f, "  : {:0>16b} {}", pencil_mask, opts_left)?,
            }
        }
        write!(f, "{} cells unpenned", self.num_unpenned_cells)
    }
}

/// Generates a clue list from a string (where `'.'` or `'0'` represent an empty cell and
/// `"123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"` is used as a case-insensitive lookup for grid labels)
pub fn clues_from_str(s: &str) -> Vec<Option<usize>> {
    s.chars()
        .filter_map(|c| match c {
            '.' | '0' => Some(None),
            '1'..='9' => Some(Some(c as usize - '1' as usize)),
            'A'..='Z' => Some(Some(c as usize - 'A' as usize + 9)),
            'a'..='z' => Some(Some(c as usize - 'a' as usize + 9)),
            _ => None,
        })
        .collect_vec()
}

/// The possible ways that a sudoku solve can fail
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum Error {
    /// The number of clues doesn't match the number of cells in the grid
    ClueLenMismatch { clue_len: usize, num_cells: usize },
    /// No solutions were found (i.e. the grid is insoluble)
    NoSolutions,
    /// Multiple solutions were found (i.e. the grid is ambiguous)
    MultipleSolutions,
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::ClueLenMismatch {
                clue_len,
                num_cells,
            } => write!(
                f,
                "Length mismatch: clues has length {}, grid has {} cells",
                clue_len, num_cells
            ),
            Error::NoSolutions => write!(f, "Puzzle has no solutions"),
            Error::MultipleSolutions => write!(f, "Puzzle is ambiguous (has multiple solutions)"),
        }
    }
}

impl std::error::Error for Error {}