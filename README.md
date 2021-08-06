# Kuudos

Creation tools for spicy sudokus.

Spicy Sudokus?  What's spicy about them?

### Non-spicy sudokus

Let's start with some non-spicy sudokus (the [Korma](https://en.wikipedia.org/wiki/Korma) of sudoku,
if you will).  Everyone knows about our venerable friend, the classic 9x9 sudoku:

\<insert 9x9 sudoku and soln\>

But most sudoku enthusiasts will know that this isn't the only possible sudoku shape.  We can also
create 4x4 and 6x6 sudokus (i.e. the boxes are 2x2 and 2x3 respectively).  However, these smaller
puzzles are generally easy to the point of being boring (there just aren't enough degrees of
freedom).  So easy, in fact, that I don't feel bad about not providing solutions :wink::

\<insert 4x4 and 6x6 sudoku\>

There are more square sudokus (16x16, 12x12, etc.), but they tend to be tedious to solve and aren't
actually very spicy (perhaps something like a [Pasanda](https://en.wikipedia.org/wiki/Pasanda)).

### Getting spicy

Here's the thing: there's nothing about the definition of a sudoku that _requires_ it to be square.
Obviously, square sudokus are easy to print and visualise, but squareness is an arbitrary limitation
imposed by humans.  But computers have no such qualms.  And since this project is about programming
computers, let's say goodbye to arbitrary restrictions and see what comes out!

At the end of the day, a sudoku is simply a set of cells.  These cells are split into _groups_
(usually rows, columns or boxes) where digits can't be placed into two different cells which share a
group.  For example, our venerable 9x9 has 81 cells and 27 groups (9 rows, 9 columns and 9 boxes).
But why limit ourselves to square grids?  Surely any set of equally sized groups is a sudoku?

Almost.  Well technically, yes - by our definition they are all sudokus.  But as humans, it would be
nice if we could actually print the sudokus (onto paper or just the screen) for us to solve.

It's worth noting at this point that Kuudos itself doesn't care what we give it - it will happily
generate and solve puzzles with any shape you feel like throwing at it - but it may be next to
impossible for a human to digest the results if they can't be printed.

### Re-arranging boxes

The easiest way to spice things up a little is to keep grouping our cells into boxes, but allow
these boxes to be arranged in something other than a square grid.  Thus, boxes now become
quadrilaterals and we have already unlocked new possibilities, such as the 'star':

\<insert 5-pointed star\>

Here, each row/column is slightly curved, but the rules are still the same: digits can't be repeated
in a box, and rows and columns are generated by following the opposite sides of edges.  More complex
and interesting shapes are also possible, such as the 'tear-drop':

\<insert teardrop\>

Whether or not these shapes are _easy_ to solve is debatable, but they are at least semi-intuitive
for humans.  I think that these are something equivalent to
[Tikka](https://en.wikipedia.org/wiki/Tikka_\(food\)) - a bit of a tingle for the little grey cells
but really nothing too spicy.

As a fun aside, it turns out that both of these shapes are incredibly over-defined.  What I mean by
this is that there are very few ways to actually fill the cells with digits, without violating one
of the group constraints.  This basically means that they are extremely easy for computers to solve,
but painfully difficult for us poor humans.  Also, the puzzle generator in Kuudos tries to make
interesting puzzles by taking away as many clues as possible whilst still having exactly one
solution.  I think you can see where this is going - the solver would merrily take out loads and
loads of clues (everything's fine, there's still one solution) and the result is nearly impossible
for any relatively inexperienced sudoku player like me.

### Splitting boxes

Anyway, we're not done yet!  In fact, we've barely got started...
