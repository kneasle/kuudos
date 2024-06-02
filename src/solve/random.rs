use itertools::Itertools;
use rand::{prelude::SliceRandom, Rng};

use crate::Shape;

use super::{
    naive::Config,
    partial::{Partial, Table},
    Error, Grid, RandomSolver, Solution, Solver,
};

/// A [`Solver`] similar to [`Naive`] that chooses new digits randomly, thus generating a random
/// solution.  This randomness is only really useful for populating empty puzzle grids - for
/// standard solvers it only causes a slow-down.
#[derive(Debug, Clone)]
pub struct NaiveRandom<'s> {
    shape: &'s Shape,
    table: Table,
}

impl<'s> NaiveRandom<'s> {
    fn recursive_solve(&self, mut partial: Partial, rng: &mut impl Rng) -> Result<Partial, Error> {
        // Solve all 'naked singles' - i.e. pen any cell which has only one possible digit.  These
        // are cheap to compute and don't require backtracking, so we perform them in-place on the
        // current Partial.
        partial.solve_naked_singles(&self.table);
        // If we already have a solution, then return with success
        if partial.is_solved() {
            return Ok(partial);
        }
        // If `partial` isn't solved but no more naked singles are found, then we pick a non-penned
        // cell and speculatively try each value to see which solutions are found.  We always want
        // to remove as much of the possible search space as possible, so we pick (one of) the
        // cells with the fewest options
        let cell_to_guess = partial.least_pencilled_cell();
        // Get and shuffle the list of available digits
        let mut digits_to_guess = partial
            .pencilled_digit_iter(cell_to_guess)
            .unwrap()
            .collect_vec();
        digits_to_guess.shuffle(rng);
        for digit in digits_to_guess {
            // Create a fresh `Partial` and guess this digit
            let mut new_partial = partial.clone();
            new_partial.pen(&self.table, cell_to_guess, digit);
            let solution_to_new_partial = self.recursive_solve(new_partial, rng);
            // Deal with the result of this speculative search
            match solution_to_new_partial {
                // If the search branch found a solution, then so did we
                Ok(solution_partial) => return Ok(solution_partial),
                // If the search branch found no solutions, then keep searching
                Err(Error::NoSolutions) => {}
                // If the iteration limit was exceeded, then terminate the search
                Err(Error::IterationLimitReached) => return Err(Error::IterationLimitReached),
                // All other error states are unreachable
                Err(Error::ClueLenMismatch { .. }) | Err(Error::MultipleSolutions) => {
                    unreachable!()
                }
            }
        }
        // If all possible digits were tried in `cell_to_guess` and none generated a solution, then
        // we have fully explored the possible ways of solving `partial` without finding a
        // solution.  Therefore, we declare that `partial` has no solution.
        Err(Error::NoSolutions)
    }
}

impl<'s> Solver<'s> for NaiveRandom<'s> {
    type Config = Config;

    fn new(shape: &'s Shape, _config: Config) -> Self {
        Self {
            shape,
            table: Table::from_shape(shape),
        }
    }
}

impl<'s> RandomSolver<'s> for NaiveRandom<'s> {
    fn solve_random(&self, clues: &Grid, rng: &mut impl Rng) -> Result<Solution, Error> {
        // Create a partial solution where only the given clues are penned
        let unsolved_partial = Partial::from_clues(&self.table, self.shape, clues)?;
        // Run recursive backtracking on this grid, and extract the solution or propagate the error
        self.recursive_solve(unsolved_partial, rng)
            .map(|solved_partial| {
                // Extract the digits, panicking if somehow the grid is still unsolved
                solved_partial
                    .into_solved_digits()
                    .expect("Solver returned unsolved grid")
            })
    }
}
