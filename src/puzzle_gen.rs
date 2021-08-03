use itertools::Itertools;
use rand::prelude::*;

use crate::{
    indexed_vec::{CellIdx, IdxType},
    solve::{MultipleSolnSolver, Solver},
    Shape, Symmetry,
};

type Grid = Vec<Option<usize>>;
type Solution = Vec<usize>;

/// A generator for random sudoku puzzles.
#[derive(Debug, Clone)]
pub struct PuzzleGen<'shp, 'symm, F: Solver<'shp>, S: MultipleSolnSolver<'shp>> {
    shape: &'shp Shape,
    symmetry: &'symm Symmetry,
    config: Config,
    grid_filler: F,
    puzzle_solver: S,
}

impl<'shp, 'symm, F: Solver<'shp>, S: MultipleSolnSolver<'shp>> PuzzleGen<'shp, 'symm, F, S> {
    /// Creates a new puzzle generator for a given [`Shape`]
    pub fn new(
        shape: &'shp Shape,
        symmetry: &'symm Symmetry,
        grid_filler_config: F::Config,
        solver_config: S::Config,
        config: Config,
    ) -> Self {
        Self {
            shape,
            symmetry,
            config,
            grid_filler: F::new(shape, grid_filler_config),
            puzzle_solver: S::new(shape, solver_config),
        }
    }

    /// Generate a grid, returning `None` if the solvers always failed
    pub fn generate(&self) -> Option<(Grid, Solution)> {
        let empty_grid: Grid = vec![None; self.shape.num_cells()];
        let mut rng = thread_rng();
        let mut best_puzzle: Option<(usize, Grid, Solution)> = None;
        // Repeatedly fill the empty grid
        for _ in 0..self.config.num_grids {
            let filled_grid = self.grid_filler.solve(&empty_grid).ok()?;
            self.gen_puzzles_from_grid(&filled_grid, &mut rng, &mut best_puzzle);
        }
        // Extract the best grid and return it if it exists
        best_puzzle.map(|(_score, grid, soln)| (grid, soln))
    }

    /// Remove clues from the grid to generate a puzzle
    fn gen_puzzles_from_grid(
        &self,
        filled_grid: &[usize],
        rng: &mut ThreadRng,
        best_puzzle: &mut Option<(usize, Grid, Solution)>,
    ) {
        // Equivalence classes nearer to the start of this list will be removed first
        let equiv_classes = {
            let mut cs = self.symmetry.equiv_classes();
            cs.shuffle(rng);
            cs
        };
        let mut clues = filled_grid.iter().copied().map(Some).collect_vec(); // Start with a full grid
        let mut num_solver_runs = 0; // Accumulator to track how many times the solver has been run

        self.recursive_remove_clues(
            filled_grid,
            &equiv_classes,
            &mut clues,
            &mut num_solver_runs,
            best_puzzle,
            0,
        );
    }

    /// Attempt to reduce the number of clues in the puzzle until the grid becomes unsolvable, with
    /// some amount of backtracking to search for other solutions
    fn recursive_remove_clues(
        &self,
        filled_grid: &[usize],
        equiv_classes: &[Vec<CellIdx>],
        clues: &mut [Option<usize>],
        num_solver_runs: &mut usize,
        best_puzzle: &mut Option<(usize, Grid, Solution)>,
        depth: usize,
    ) {
        // Early return if there are no more equivalence classes to remove (i.e. we've reached the
        // max recursion depth)
        if depth == filled_grid.len() {
            return;
        }
        // Stop the search if we've run the solver too many times
        if *num_solver_runs >= self.config.num_solver_runs_per_grid {
            return;
        }
        // Remove the clues contained in the `depth`th equivalence class
        for cell_idx in &equiv_classes[depth] {
            clues[cell_idx.to_idx()] = None;
        }
        // Test whether or not the grid is soluble
        *num_solver_runs += 1;
        if self.puzzle_solver.solve_multiple_solns(clues).is_ok() {
            // Update the `best_puzzle` if this is better
            let score = clues.iter().filter(|v| v.is_some()).count();
            match best_puzzle {
                Some((best_score, best_grid, best_soln)) => {
                    if score < *best_score {
                        // If this better than the best grid/solution, then this becomes the best
                        // grid/solution
                        best_grid.clear();
                        best_grid.extend_from_slice(clues);
                        best_soln.clear();
                        best_soln.extend_from_slice(filled_grid);
                        *best_score = score;
                    }
                }
                // If no grids have been discovered, then this one is automatically the best
                None => *best_puzzle = Some((score, clues.to_vec(), filled_grid.to_vec())),
            }
            // If this grid _is_ soluble, then try taking away more equivalence classes
            self.recursive_remove_clues(
                filled_grid,
                equiv_classes,
                clues,
                num_solver_runs,
                best_puzzle,
                depth + 1,
            );
        }
        // Re-add the clues we removed
        for cell_idx in &equiv_classes[depth] {
            clues[cell_idx.to_idx()] = Some(filled_grid[cell_idx.to_idx()]);
        }
        // Now that we've added these cells back, we try to remove the next equivalence class
        self.recursive_remove_clues(
            filled_grid,
            equiv_classes,
            clues,
            num_solver_runs,
            best_puzzle,
            depth + 1,
        );
    }
}

/// Configuration parameters for a [`PuzzleGen`]erator
#[derive(Debug, Clone)]
pub struct Config {
    /// How many solved grids are generated
    num_grids: usize,
    /// For each solved grid, how many candidate puzzles are generated
    num_solver_runs_per_grid: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            num_grids: 10,
            num_solver_runs_per_grid: 200,
        }
    }
}
