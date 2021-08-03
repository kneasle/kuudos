#![allow(dead_code)]

use itertools::Itertools;
use kuudos::{
    puzzle_gen::{self, PuzzleGen},
    shape::RenderingOpts,
    solve::{
        clues_from_str,
        naive::{self, Naive},
        Solver,
    },
    Shape,
};

const VALUE_NAMES: &str = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

fn main() {
    let (shape, clues, soln) = if true {
        // let s = Shape::star2x2(5, PI);
        let builder = shapes::nine_mens_morris().unwrap();
        std::fs::write("bdr.svg", builder.as_svg(40.0)).unwrap();

        let (s, symm) = builder.build().unwrap();
        let (clues, soln) = PuzzleGen::<Naive, Naive>::new(
            &s,
            &symm,
            naive::Config::default(),
            naive::Config::default(),
            puzzle_gen::Config::default(),
        )
        .generate()
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
            .solve(&clues)
            .unwrap();
        (s, clues, soln)
    };

    // Write SVG files of the shape's clues, cell numbering and solution
    let clue_vec = clues
        .iter()
        .map(|x| x.map(|v| VALUE_NAMES.chars().nth(v).unwrap()))
        .collect_vec();
    let soln_vec = soln
        .iter()
        .map(|digit| VALUE_NAMES.chars().nth(*digit).unwrap())
        .map(Some)
        .collect_vec();
    let cell_name_vec = VALUE_NAMES
        .chars()
        .cycle()
        .take(shape.num_cells())
        .map(Some)
        .collect_vec();
    // Generate all the files
    let write_svg_str = |clues: Vec<Option<char>>, path: &str| {
        let svg_str = shape.svg_string(&RenderingOpts::default(), 40.0, &clues);
        std::fs::write(path, svg_str).unwrap();
    };
    write_svg_str(clue_vec, "puzzle.svg");
    write_svg_str(soln_vec, "soln.svg");
    write_svg_str(cell_name_vec, "cells.svg");
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum DisplayType {
    Clues,
    Solution,
    CellNames,
}

mod shapes {
    use angle::Deg;
    use kuudos::{
        builder::{BoxAddError, Builder, EdgeLinkStyle},
        regular_polygon_inradius, Side, V2Ext, V2,
    };

    pub fn nine_mens_morris() -> Result<Builder, BoxAddError> {
        const LINK_DEPTH: f32 = 0.2;

        let mut bdr = Builder::new(3, 3, 4);

        // Add line(s) of connected cells
        let line_box_1 = bdr.add_box_square(V2::UP * 3.0, 1.0, Deg(0.0));
        let line_box_2 = bdr.extrude_edge_with_opts(line_box_1, Side::Top, 1.0, LINK_DEPTH)?;
        let line_box_3 = bdr.extrude_edge_with_opts(line_box_2, Side::Top, 1.0, LINK_DEPTH)?;
        // Add diagonal line(s) of 'floating' cells
        let up_right = V2::UP + V2::RIGHT;
        bdr.add_box_square(up_right * 3.0, 1.0, Deg(0.0));
        let corner_box_2 = bdr.add_box_square(up_right * (6.0 + LINK_DEPTH), 1.0, Deg(0.0));
        let corner_box_3 = bdr.add_box_square(up_right * (9.0 + LINK_DEPTH * 2.0), 1.0, Deg(0.0));
        // Link the lines of boxes
        bdr.link_edges(
            line_box_2,
            Side::Right,
            corner_box_2,
            Side::Left,
            EdgeLinkStyle::Linear,
        );
        bdr.link_edges(
            line_box_3,
            Side::Right,
            corner_box_3,
            Side::Left,
            EdgeLinkStyle::Linear,
        );
        bdr.link_edges(
            corner_box_2,
            Side::Bottom,
            bdr.rotational_copy_of(line_box_2, 1).unwrap(),
            Side::Left,
            EdgeLinkStyle::Linear,
        );
        bdr.link_edges(
            corner_box_3,
            Side::Bottom,
            bdr.rotational_copy_of(line_box_3, 1).unwrap(),
            Side::Left,
            EdgeLinkStyle::Linear,
        );

        Ok(bdr)
    }

    pub fn race_track() -> Result<Builder, BoxAddError> {
        // Create a Builder with 3-fold symmetry
        let mut bdr = Builder::new(3, 3, 3);

        // Create a new parallelogram box to make a hexagon
        let hex_box =
            bdr.add_box_parallelogram(V2::ZERO, V2::UP, bdr.rotate_point_by_steps(V2::UP, 1))?;
        // Extrude a square face off one side
        let outer_box = bdr.extrude_edge(hex_box, Side::Top)?;
        // Link the outer boxes together
        bdr.link_edges(
            outer_box,
            Side::Right,
            bdr.rotational_copy_of(outer_box, 1).unwrap(),
            Side::Left,
            EdgeLinkStyle::Arc,
        );

        // Rotate the hexagon so that the side we extruded faces upwards
        bdr.rotate(Deg(-30.0));
        Ok(bdr)
    }

    pub fn radio_waves() -> Result<Builder, BoxAddError> {
        let triangle_radius = regular_polygon_inradius(3, 3.0);

        let mut bdr = Builder::new(3, 3, 3);
        // Inner ring
        let bottom_box = bdr.add_box_square(V2::UP * (triangle_radius + 1.5), 1.0, Deg(0.0));
        let _corner_box = bdr.add_box_parallelogram(
            V2::RIGHT * 1.5 + V2::UP * triangle_radius,
            V2::UP,
            V2::UP.rotate(Deg(120.0)),
        )?;
        // Outer rings
        let stack_box_1 = bdr.extrude_edge(bottom_box, Side::Top)?;
        let stack_box_2 = bdr.extrude_edge(stack_box_1, Side::Top)?;
        for stack_box in [stack_box_1, stack_box_2] {
            bdr.link_edges(
                stack_box,
                Side::Right,
                bdr.rotational_copy_of(stack_box, 1).unwrap(),
                Side::Left,
                EdgeLinkStyle::Arc,
            )
            .unwrap();
        }

        Ok(bdr)
    }

    pub fn cube() -> Result<Builder, BoxAddError> {
        /* Cube net layout and box naming:
         *
         *              top_box
         *                 |
         * left_box -- front_box -- right_box
         *                 |
         *            bottom_box
         *                 |
         *             back_box
         */
        let mut bdr = Builder::new(4, 4, 1);
        // Add boxes
        let front_box = bdr.add_box_square(V2::ZERO, 1.0, Deg(0.0));
        let top_box = bdr.extrude_edge(front_box, Side::Top)?;
        let left_box = bdr.extrude_edge(front_box, Side::Left)?;
        let right_box = bdr.extrude_edge(front_box, Side::Right)?;
        let bottom_box = bdr.extrude_edge(front_box, Side::Bottom)?;
        let back_box = bdr.extrude_edge(bottom_box, Side::Bottom)?;

        // Link the top-right-left-bottom cycle
        bdr.link_edges(
            top_box,
            Side::Right,
            right_box,
            Side::Top,
            EdgeLinkStyle::Hidden,
        )
        .unwrap();
        bdr.link_edges(
            top_box,
            Side::Left,
            left_box,
            Side::Top,
            EdgeLinkStyle::Hidden,
        )
        .unwrap();
        bdr.link_edges(
            left_box,
            Side::Bottom,
            bottom_box,
            Side::Left,
            EdgeLinkStyle::Hidden,
        )
        .unwrap();
        // Link the front-right-back-left cycle
        bdr.link_edges(
            right_box,
            Side::Right,
            back_box,
            Side::Right,
            EdgeLinkStyle::Hidden,
        )
        .unwrap();

        Ok(bdr)
    }

    pub fn useless_wheel() -> Result<Builder, BoxAddError> {
        const LINK_DEPTH: f32 = 1.0;

        let mut bdr = Builder::new(3, 3, 4);
        // We build the top-right corner, and let the symmetry expand the rest
        let central_box = bdr.add_box_parallelogram(V2::ZERO, V2::UP, V2::RIGHT)?;
        let right_box = bdr.extrude_edge_with_opts(central_box, Side::Right, 1.0, LINK_DEPTH)?;
        let top_box = bdr.extrude_edge_with_opts(right_box, Side::Top, 1.0, LINK_DEPTH)?;
        bdr.link_edges(
            top_box,
            Side::Left,
            bdr.rotational_copy_of(right_box, -1).unwrap(),
            Side::Bottom,
            EdgeLinkStyle::Linear,
        )
        .unwrap();

        Ok(bdr)
    }

    pub fn diamonds() -> Result<Builder, BoxAddError> {
        let mut bdr = Builder::new(4, 4, 4);
        // Add 2x2 cells in the middle, and diamonds in the corners
        let tr_central_box = bdr.add_box_parallelogram(V2::ZERO, V2::UP, V2::RIGHT)?;
        let right_diamond = bdr.add_box_square_by_corner(V2::RIGHT * 4.0, 1.0, Deg(45.0));
        // Link edges
        bdr.link_edges(
            tr_central_box,
            Side::Right,
            right_diamond,
            Side::Left,
            EdgeLinkStyle::Arc,
        )
        .unwrap();
        // Link edges
        bdr.link_edges(
            right_diamond,
            Side::Bottom,
            bdr.rotational_copy_of(tr_central_box, 1).unwrap(),
            Side::Top,
            EdgeLinkStyle::Arc,
        )
        .unwrap();

        Ok(bdr)
    }

    pub fn mini_4x4() -> Result<Builder, BoxAddError> {
        let mut bdr = Builder::new(4, 4, 4);
        // Add 2x2 cells in the middle (again, we're building the top-right corner)
        let central_box = bdr.add_box_parallelogram(V2::ZERO, V2::UP, V2::RIGHT)?;
        bdr.link_edges(
            central_box,
            Side::Right,
            bdr.rotational_copy_of(central_box, 1).unwrap(),
            Side::Top,
            EdgeLinkStyle::Arc,
        )
        .unwrap();

        Ok(bdr)
    }

    pub fn space_station() -> Result<Builder, BoxAddError> {
        let inradius_of_triangle = regular_polygon_inradius(3, 2.0);

        let mut bdr = Builder::new(3, 2, 3);
        // Add 2x2 cells in the middle (again, we're building the top-right corner)
        let horiz_box =
            bdr.add_box_square_by_corner(V2::new(inradius_of_triangle, 1.0), 1.0, Deg(0.0));
        let _upper_box = bdr.add_box_parallelogram(
            V2::new(inradius_of_triangle, -1.0),
            V2::UP.rotate(Deg(30.0)),
            V2::RIGHT,
        )?;
        let _lower_box = bdr.add_box_parallelogram(
            V2::new(inradius_of_triangle, 1.0),
            V2::DOWN.rotate(Deg(-30.0)),
            V2::RIGHT,
        )?;

        let orbiting_face = bdr.extrude_edge_with_opts(horiz_box, Side::Right, 1.0, 1.4)?;
        bdr.link_edges(
            orbiting_face,
            Side::Bottom,
            bdr.rotational_copy_of(orbiting_face, 1).unwrap(),
            Side::Top,
            EdgeLinkStyle::Arc,
        )
        .unwrap();

        Ok(bdr)
    }

    pub fn rot_sym_classic() -> Result<Builder, BoxAddError> {
        let mut bdr = Builder::new(3, 3, 4);
        let centre_box = bdr.add_box_square(V2::ZERO, 1.0, Deg(0.0));
        let right_box = bdr.extrude_edge(centre_box, Side::Right)?;
        let _top_box = bdr.extrude_edge(right_box, Side::Top)?;
        Ok(bdr)
    }
}
