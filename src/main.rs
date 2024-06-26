use kuudos::{
    image::RenderingOpts,
    puzzle_gen::{self, PuzzleGen},
    shape::{examples, CellVec, Shape},
    solve::{
        clues_from_str,
        naive::{self, Naive},
        random::NaiveRandom,
        SingleSolnSolver, Solver,
    },
    Symmetry,
};

const VALUE_NAMES: &str = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

fn main() {
    let (shape, clues, soln) = if true {
        // let s = Shape::star2x2(5, PI);
        let (s, symm) = if true {
            examples::star(5)
        } else {
            let shape = Shape::square(3, 2);
            let symm = Symmetry::asymmetric(&shape);
            (shape, symm)
        };

        let (clues, soln) = PuzzleGen::<NaiveRandom, Naive>::new(
            &s,
            &symm,
            naive::Config::default(),
            naive::Config::default(),
            puzzle_gen::Config::default(),
        )
        .generate_from_seed(0)
        .unwrap();
        (s, clues, soln)
    } else {
        let s = Shape::classic();
        let clues = clues_from_str(
            // "004300209 005009001 070060043
            //  006002087 190007400 050083000
            //  600000105 003508690 042910300",
            // "..24..... .......98 ..8..16..
            //  .5....826 ....4.7.. 7..8...5.
            //  ..167.... 3...92... ..9......",
            "4.....8.5 .3....... ...7.....
             .2.....6. ....8.4.. ....1....
             ...6.3.7. 5..2..... 1.4......",
        );
        let soln = Naive::new(&s, naive::Config::default())
            .solve_single_soln(&clues)
            .unwrap();
        (s, clues, soln)
    };

    // Write SVG files of the shape's clues, cell numbering and solution
    let clue_vec = clues
        .iter()
        .map(|x| x.map(|v| VALUE_NAMES.chars().nth(v).unwrap()))
        .collect::<CellVec<_>>();
    let soln_vec = soln
        .iter()
        .map(|digit| VALUE_NAMES.chars().nth(*digit).unwrap())
        .map(Some)
        .collect::<CellVec<_>>();
    let cell_name_vec = VALUE_NAMES
        .chars()
        .cycle()
        .take(shape.num_cells())
        .map(Some)
        .collect::<CellVec<_>>();
    // Generate all the files
    let write_svg_str = |clues: CellVec<Option<char>>, path: &str| {
        let svg_str = shape.svg_string(&RenderingOpts::default(), 40.0, &clues);
        std::fs::write(path, svg_str).unwrap();
    };
    write_svg_str(clue_vec, "puzzle.svg");
    write_svg_str(soln_vec, "soln.svg");
    write_svg_str(cell_name_vec, "cells.svg");
}
