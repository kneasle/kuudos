use std::fmt::{Display, Formatter};

use itertools::Itertools;

use crate::Shape;

pub mod naive;

pub trait Solver<'s> {
    /// Create a new solver that can be used to solve grids
    fn new(shape: &'s Shape) -> Self;

    /// Solve a puzzle
    fn solve(&self, clues: &[Option<usize>]) -> Result<Vec<usize>, Error>;
}

pub trait SolverMultipleSolns<'s>: Solver<'s> {
    /// Solve a puzzle, checking for multiple solutions
    fn solve_multiple_solns(&self, clues: &[Option<usize>]) -> Result<Vec<usize>, Error>;
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
