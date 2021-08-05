use crate::Shape;

use super::{
    partial::{Partial, Table},
    Error, MultipleSolnSolver, Solver, WithDifficulty,
};

/// A pretty naive solver, using only 'naked singles' and backtracking as well as relying on
/// runtime branches for all configuration.  Nevertheless, this naive approach is fast enough to
/// efficiently solve most human-friendly puzzles.  It is, however, extremely bad at detecting
/// unsolvable grids, making it unsuitable for tasks such as populating partially complete grids.
/// It is mainly intended as a cheap-and-cheerful test rig for other parts of the code.
#[derive(Debug, Clone)]
pub struct Naive<'s> {
    shape: &'s Shape,
    config: Config,
    table: Table,
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
        let partial = Partial::from_clues(&self.table, self.shape, clues)?;
        // Run recursive backtracking on this grid, and extract the solution or propagate the error
        self.recursive_solve(partial, check_multiple_solns, &mut num_digits_penned)
            // Extract the solution
            .map(|partial| {
                partial
                    .into_solved_digits()
                    .expect("Solver returned unsolved grid")
            })
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
        let num_naked_singles_solved = partial.solve_naked_singles(&self.table);
        *num_digits_penned += num_naked_singles_solved;

        // If the grid has not been solved by naked singles, then it is either solved or a branch
        // in the search tree is required

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
        let cell_to_guess = partial.least_pencilled_cell();

        // Place to store one solution.  If we are checking for multiple solutions, the solver will
        // place the first solution in here, then search for a second one.  If a second one is
        // later found, then `solution` will be `Some(p)` and the search can be aborted.
        let mut solution: Option<Partial> = None;
        let mut unexplored_pencil_mask = partial.get_pencil_mask(cell_to_guess).unwrap();
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
            new_partial.pen(&self.table, cell_to_guess, assumed_digit);
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
        Self {
            shape,
            config,
            table: Table::from_shape(shape),
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
