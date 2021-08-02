use std::fmt::{Display, Formatter};

use itertools::Itertools;

use crate::{
    indexed_vec::{CellIdx, CellVec, IdxType},
    Shape,
};

use super::{Error, MultipleSolnSolver, Solver, WithDifficulty};

/// A pretty naive solver, using only 'naked singles' and backtracking.  Nevertheless, this naive
/// approach is fast enough to efficiently solve most human-friendly puzzles.  It is, however,
/// extremely bad at detecting unsolvable grids, making it unsuitable for tasks such as populating
/// partially complete grids.  It is mainly intended as a cheap-and-cheerful test rig for other
/// parts of the code.
#[derive(Debug, Clone)]
pub struct Naive<'s> {
    shape: &'s Shape,
    config: Config,
    /// For a cell with index `i`, `affected_cells[i]` contains the set of cells which share a
    /// house with cell `i` (**excluding** `i` itself)
    affected_cells: CellVec<Vec<CellIdx>>,
}

impl<'s> Naive<'s> {
    fn solve_(
        &self,
        clues: &[Option<usize>],
        check_multiple_solns: bool,
    ) -> Result<(Vec<usize>, f32), Error> {
        if self.shape.num_cells() != clues.len() {
            return Err(Error::ClueLenMismatch {
                clue_len: clues.len(),
                num_cells: self.shape.num_cells(),
            });
        }

        let mut num_digits_penned = 0; // We use number of cells filled as our metric for
                                       // difficulty

        // Create a partial solution where only the given clues are penned
        let mut partial = Partial::empty(self.shape);
        clues
            .iter()
            .enumerate()
            .map(|(idx, v)| (CellIdx::from_idx(idx), v))
            .filter_map(|(i, maybe_v)| maybe_v.map(|v| (i, v)))
            .for_each(|(cell_idx, value)| partial.pen(self, cell_idx, value));

        // Run recursive backtracking on this grid, and extract the solution or propagate the error
        self.recursive_solve(partial, check_multiple_solns, &mut num_digits_penned)
            .map(Partial::into_solved_digits) // Extract the solution
            .map(|grid| (grid, num_digits_penned as f32)) // Add the difficulty to the return value
    }

    fn recursive_solve(
        &self,
        mut partial: Partial,
        check_multiple_solns: bool,
        num_digits_penned: &mut usize,
    ) -> Result<Partial, Error> {
        // Repeatedly pen in 'naked singles' - i.e. a cell with only one possible digit.  These are
        // cheap to compute and don't require backtracking, so we perform them in-place on the
        // current Partial.
        loop {
            let mut has_naked_singles = false;
            for (cell_idx, _cell) in self.affected_cells.indexed_iter() {
                // If this cell is a naked single then we pen its value
                if let Some(only_digit) = partial.naked_single_digit(cell_idx) {
                    has_naked_singles = true;
                    partial.pen(self, cell_idx, only_digit);
                    *num_digits_penned += 1;
                    // Whenever we pen a cell, check that the invariants are still upheld (but only
                    // in debug mode)
                    partial.debug_assert_invariants();
                }
            }
            // Exit the loop when no more naked singles are possible
            if !has_naked_singles {
                break;
            }
        }

        // If the grid has not been solved by naked singles, then it is either solved or a branch
        // is required

        // If the value of `partial` passed to this function is soluble using only naked singles
        // then it permits precisely one solution, which `partial` now contains.  So we return that
        // single solution
        if partial.is_solved() {
            return Ok(partial);
        }

        // If `partial` isn't solved but no more naked singles are found, then we pick a non-penned
        // cell and speculatively try each value to see which solutions are found.  We always want
        // to remove as much of the possible search space as possible, so we pick (one of) the
        // cells with the fewest options
        let min_idx = partial
            .num_pencil_marks
            .indexed_iter()
            .min_by_key(|(_cell_idx, num_opts)| **num_opts)
            .map(|(cell_idx, _)| cell_idx)
            .expect("Grid has no cells");

        let mut solution: Option<Partial> = None;
        let mut unexplored_pencil_mask = partial.pencil_masks[min_idx];
        while unexplored_pencil_mask != 0 {
            // Repeatedly pick the highest pencil mark in the cell, and remove its flag bit
            let assumed_digit = 63 - unexplored_pencil_mask.leading_zeros() as usize;
            unexplored_pencil_mask &= !(1 << assumed_digit);
            // Add +1 difficulty because we had to guess.  This is a pretty terrible difficulty
            // metric - the `solve::human::Human` solver implements a much better approximation of
            // human difficulty
            *num_digits_penned += 1;
            if *num_digits_penned >= self.config.max_iterations {
                return Err(Error::IterationLimitReached);
            }
            // Create a new `Partial` which assumes that this digit is taken
            let mut new_partial = partial.clone();
            new_partial.pen(self, min_idx, assumed_digit);
            match self.recursive_solve(new_partial, check_multiple_solns, num_digits_penned) {
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
                // If no solutions where found, then we keep looking for a solution
                Err(Error::NoSolutions) => {}
                // If the solver produced too many solutions, then this whole branch is unsolvable
                // (because the other branches can only generate more solutions)
                Err(Error::MultipleSolutions) => return Err(Error::MultipleSolutions),
                // Reaching the iteration limit forces the search to terminate, and therefore
                // should make every recursive call return
                Err(Error::IterationLimitReached) => return Err(Error::IterationLimitReached),

                // This error can only be generated by the top-level `solve()` function
                Err(Error::ClueLenMismatch { .. }) => unreachable!(),
            }
        }

        // If multiple solutions exist down this branch, we would have taken an early return.
        // Therefore, this branch contains either 0 or 1 solutions (encoded by `solution` being
        // `None` or `Some(_)` respectively).
        solution.ok_or(Error::NoSolutions)
    }
}

impl<'s> Solver<'s> for Naive<'s> {
    type Config = Config;

    fn new(shape: &'s Shape, config: Config) -> Self {
        let mut affected_cells: CellVec<Vec<CellIdx>> =
            CellVec::repeat(Vec::new(), shape.num_cells());
        for group in &shape.groups {
            for &cell_idx in &group.cells {
                affected_cells[cell_idx].extend(group.cells.iter().filter(|idx| **idx != cell_idx));
            }
        }
        // Deduplicate the lists of affected cells
        for cells in affected_cells.iter_mut() {
            cells.sort_unstable();
            cells.dedup();
        }
        Self {
            shape,
            config,
            affected_cells,
        }
    }

    fn solve(&self, clues: &[Option<usize>]) -> Result<Vec<usize>, Error> {
        self.solve_(clues, false) // Delegate to the generic solver
            .map(|(soln, _difficulty)| soln) // Discard the difficulty
    }
}

impl<'s> MultipleSolnSolver<'s> for Naive<'s> {
    fn solve_multiple_solns(&self, clues: &[Option<usize>]) -> Result<Vec<usize>, Error> {
        self.solve_(clues, true) // Delegate to the generic solver
            .map(|(soln, _difficulty)| soln) // Discard the difficulty
    }
}

impl<'s> WithDifficulty<'s> for Naive<'s> {
    fn solve_with_difficulty(&self, clues: &[Option<usize>]) -> Result<(Vec<usize>, f32), Error> {
        self.solve_(clues, true) // Delegate to the generic solver
    }
}

/// Configuration parameters for the [`Naive`] solver
#[derive(Debug, Clone)]
pub struct Config {
    /// Number of iterations (i.e. branches) that the [`Naive`] solver will take before rejecting
    /// the grid as unsolvable
    max_iterations: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            max_iterations: 10_000_000, // 10M iterations sounds like enough
        }
    }
}

/// The state of a `Partial`ly solved sudoku grid
#[derive(Debug, Clone)]
struct Partial {
    /// Which cells have 'penned' values
    penned_cells: CellVec<Option<usize>>,
    /// How many cells are 'unpenned'.  A `Partial` is solved if this is zero.
    ///
    /// Invariant: `num_unpenned_cells = penned_cells.iter().filter(|c| c.is_none()).count()`
    num_unpenned_cells: usize,

    /// A bitmask for each cell, where a `1` at position `i` means that digit `i` can be placed in
    /// that cell
    pencil_masks: CellVec<usize>,
    /// How many options are left per cell.
    ///
    /// Invariant: for all i: `num_options_left[i] = options_per_cell[i].count_ones()`
    num_pencil_marks: CellVec<usize>,
}

impl Partial {
    /// A `Partial` representing a completely empty grid
    fn empty(shape: &Shape) -> Partial {
        // `all_options` has `shape.num_symbols` 1s in its lowest bits, and 0s elsewhere
        let all_options = (1 << shape.num_symbols) - 1;
        Self {
            penned_cells: CellVec::repeat(None, shape.num_cells()),
            num_unpenned_cells: shape.num_cells(),
            pencil_masks: CellVec::repeat(all_options, shape.num_cells()),
            num_pencil_marks: CellVec::repeat(shape.num_symbols, shape.num_cells()),
        }
    }

    /// Pen a value in a given cell
    fn pen(&mut self, solver: &Naive, cell_idx: CellIdx, value: usize) {
        // Record the penned value in the current cell
        if self.penned_cells[cell_idx].is_none() {
            self.num_unpenned_cells -= 1;
        }
        self.penned_cells[cell_idx] = Some(value);
        // Set the number of pencil marks in this cell to be as large as possible (so this cell
        // will never get chosen for either a hidden single or as the best branch).
        self.num_pencil_marks[cell_idx] = usize::MAX;
        // 'Unpencil' this value from any cell which shares a house with this one
        for &affected_cell_idx in &solver.affected_cells[cell_idx] {
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
    fn naked_single_digit(&self, cell_idx: CellIdx) -> Option<usize> {
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
            .iter()
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
        for ((cell_idx, pencil_mask), count) in self
            .pencil_masks
            .indexed_iter()
            .zip_eq(self.num_pencil_marks.iter())
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
        for (pen, (pencil_mask, opts_left)) in self.penned_cells.iter().zip_eq(
            self.pencil_masks
                .iter()
                .zip_eq(self.num_pencil_marks.iter()),
        ) {
            match pen {
                Some(p) => writeln!(f, "{:>2}: {:0>16b} {}", p, pencil_mask, opts_left)?,
                None => writeln!(f, "  : {:0>16b} {}", pencil_mask, opts_left)?,
            }
        }
        write!(f, "{} cells unpenned", self.num_unpenned_cells)
    }
}
