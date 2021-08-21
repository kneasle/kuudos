use std::{
    collections::HashSet,
    convert::TryFrom,
    fmt::{Display, Formatter},
};

use enum_iterator::IntoEnumIterator;
use kuudos::{
    image::{self, RenderingOpts},
    indexed_vec::CellVec,
    puzzle_gen::{self, PuzzleGen},
    shape::{char_in_cell_center, examples},
    solve::{
        naive::{self, Naive},
        random::NaiveRandom,
    },
    Shape, Symmetry, V2,
};
use num_enum::{IntoPrimitive, TryFromPrimitive};
use wasm_bindgen::prelude::*;

type Digit = usize; // TODO: Make a proper `Digit` type in kuudos

const DIGIT_NAMES: [char; 35] = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
];

/// A partially solved puzzle
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct Puzzle {
    shape: Shape,
    cell_contents: CellVec<CellContents>,
    solution: CellVec<Digit>,
}

#[wasm_bindgen]
impl Puzzle {
    pub fn from_shape_id(id: usize, seed: u64) -> Option<Puzzle> {
        Self::new(PuzzleShape::from_id(id)?, seed)
    }

    fn new(puzzle_shape: PuzzleShape, seed: u64) -> Option<Puzzle> {
        let (shape, symmetry) = puzzle_shape.generate_shape();
        let puzzle_generator = PuzzleGen::<NaiveRandom, Naive>::new(
            &shape,
            &symmetry,
            naive::Config::default(),
            naive::Config::default(),
            puzzle_gen::Config::default(),
        );
        let (clues, solution) = puzzle_generator.generate_from_seed(seed)?;
        Some(Self {
            shape,
            cell_contents: clues
                .iter()
                .map(|digit| match digit {
                    Some(d) => CellContents::Clue(*d),
                    // Empty cells are pencilled with no digits
                    None => CellContents::Pencilled(HashSet::new()),
                })
                .collect(),
            solution: solution.into(),
        })
    }

    pub fn svg_string(&self, scale: f32) -> String {
        self.shape
            .svg_string_with_contents(&RenderingOpts::default(), scale, |idx, verts| {
                self.cell_contents[idx].img_elems(verts)
            })
    }
}

/// The possible contents that a cell can have whilst solving.  An empty cell is represented as
/// [`Pencilled`](CellContents::Pencilled) with an empty [`HashSet`] (i.e. no pencilled digits).
#[derive(Debug, Clone)]
enum CellContents {
    Clue(Digit),
    Pen(Digit),
    Pencilled(HashSet<Digit>),
}

impl CellContents {
    fn img_elems(&self, verts: &[V2]) -> Vec<image::Elem> {
        match self {
            CellContents::Clue(d) => {
                char_in_cell_center(image::TextStyle::Clue, verts, DIGIT_NAMES[*d])
            }
            CellContents::Pen(d) => {
                char_in_cell_center(image::TextStyle::PennedDigit, verts, DIGIT_NAMES[*d])
            }
            CellContents::Pencilled(_pencilled_digits) => todo!(),
        }
    }
}

#[derive(
    IntoEnumIterator, IntoPrimitive, TryFromPrimitive, Debug, Clone, Copy, Eq, PartialEq, Hash,
)]
#[repr(usize)]
pub enum PuzzleShape {
    Square4x4 = 0,
    Square6x6,
    Square9x9,
    Hex,
    Star,
    Triangle,
    SpaceStation,
    RaceTrack,
    Wheel,
    Mini4x4,
    // NineMensMorris,
    // RadioWaves,
    Cube,
    Diamonds,
}

impl PuzzleShape {
    fn from_id(id: usize) -> Option<Self> {
        Self::try_from(id).ok()
    }

    #[allow(dead_code)] // Not used, but included for completeness
    fn id(self) -> usize {
        self.into()
    }

    /// Generates a [`Shape`]/[`Symmetry`] pair from a `PuzzleShape`
    fn generate_shape(self) -> (Shape, Symmetry) {
        fn asymetric_square(box_width: usize, box_height: usize) -> (Shape, Symmetry) {
            let shape = Shape::square(box_width, box_height);
            let symmetry = Symmetry::asymmetric(&shape);
            (shape, symmetry)
        }

        match self {
            PuzzleShape::Square4x4 => asymetric_square(2, 2), // TODO: Pick the symmetry at random
            PuzzleShape::Square6x6 => asymetric_square(3, 2),
            PuzzleShape::Square9x9 => examples::rot_sym_classic(),
            PuzzleShape::Hex => examples::star(3), // A hexagon is a 3-pointed star
            PuzzleShape::Star => examples::star(5),
            PuzzleShape::Triangle => examples::triangle(),
            PuzzleShape::SpaceStation => examples::space_station(),
            PuzzleShape::RaceTrack => examples::race_track(),
            PuzzleShape::Wheel => examples::useless_wheel(),
            PuzzleShape::Mini4x4 => examples::mini_4x4(),
            PuzzleShape::Cube => examples::cube(),
            PuzzleShape::Diamonds => examples::diamonds(),
        }
    }
}

impl Display for PuzzleShape {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            PuzzleShape::Square4x4 => write!(f, "4x4 Square"),
            PuzzleShape::Square6x6 => write!(f, "6x6 Square"),
            PuzzleShape::Square9x9 => write!(f, "9x9 Square"),
            PuzzleShape::Hex => write!(f, "Hexagon"),
            PuzzleShape::Star => write!(f, "5-Point Star"),
            PuzzleShape::Triangle => write!(f, "Triangle"),
            PuzzleShape::SpaceStation => write!(f, "Space Station"),
            PuzzleShape::RaceTrack => write!(f, "Race Track"),
            PuzzleShape::Wheel => write!(f, "Wheel"),
            PuzzleShape::Mini4x4 => write!(f, "Mini 4x4"),
            PuzzleShape::Cube => write!(f, "Cube"),
            PuzzleShape::Diamonds => write!(f, "Diamonds"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Puzzle;

    #[test]
    fn gen_puzzle() {
        let _ = Puzzle::from_shape_id(0, 0).unwrap().svg_string(40.0);
    }
}
