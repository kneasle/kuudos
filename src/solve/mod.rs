use std::fmt::{Display, Formatter};

use itertools::Itertools;
use rand::Rng;

use crate::Shape;

pub mod naive;
mod partial;
pub mod random;

/// Type alias for a partially filled sudoku grid
pub type Grid = Vec<Option<usize>>;
/// Type alias for a fully filled (i.e. solved) sudoku grid
pub type Solution = Vec<usize>;

/// Trait implemented by all solving algorithms
pub trait Solver<'s> {
    type Config: Default;

    /// Create a new solver that can be used to solve as many grids as needed.  This function
    /// is essentially used to give the solvers an opportunity to build lookup tables for each
    /// [`Shape`].
    fn new(shape: &'s Shape, config: Self::Config) -> Self;
}

/// Trait implemented by all solvers which can check for single solutions
pub trait SingleSolnSolver<'s>: Solver<'s> {
    /// Find a single solution for a puzzle, without checking for the existence of other solutions.
    fn solve_single_soln(&self, clues: &[Option<usize>]) -> Result<Solution, Error>;
}

/// Trait implemented by all solving algorithms which check for multiple solutions
pub trait MultipleSolnSolver<'s>: Solver<'s> {
    /// Solve a puzzle, checking for multiple solutions
    fn solve(&self, clues: &[Option<usize>]) -> Result<Solution, Error>;
}

/// Trait implemented by all solving algorithms which check for multiple solutions and difficulty
/// of the solution
pub trait WithDifficulty<'s>: Solver<'s> {
    /// Solve a puzzle, checking for multiple solutions
    fn solve_with_difficulty(&self, clues: &[Option<usize>]) -> Result<(Solution, f32), Error>;
}

/// Trait implemented by all solving algorithms which check for multiple solutions and difficulty
/// of the solution
pub trait RandomSolver<'s>: Solver<'s> {
    /// Solve a puzzle for one solution, filling the cells randomly.  This is only really used for
    /// populating puzzle grids; for solving puzzles, randomness has no benefit and only slows the
    /// solver down.
    fn solve_random(&self, clues: &[Option<usize>], rng: &mut impl Rng) -> Result<Solution, Error>;
}

/// Generates a clue list from a string (where `'.'` or `'0'` represent an empty cell and
/// `"123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"` is used as a case-insensitive lookup for grid labels)
pub fn clues_from_str(s: &str) -> Grid {
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
    /// The solver hit its iteration limit before the search terminated naturally
    IterationLimitReached,
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
            Error::IterationLimitReached => write!(f, "Solver's iteration limit was reached"),
        }
    }
}

impl std::error::Error for Error {}
