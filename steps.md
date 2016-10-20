1. Open Model
2. Check links and geometry
3. Apply Loads (Std, Dead, HL-93 12ft lane IF 33% noncentrifugal load path, maximize lanes)
4. Response variables: (Moment) FXX (max or min) on flange, (Shear) FZZ on web (max AND min)
5. Run LIA
6. Load Influence Combinations: choose max/min locations to place HL-93
7. Run LSA & edit Vmax accordingly (go back to step 6)
8. Pull moment data using matlab file (top flange, web, bottom flanges)
9. Pull shear data manually using s7 (peak tool)

plate nr, location
118  A1 int
128  A1 ext
23519 P1 int
25509 P1 ext
23619 A2 int
25609 A2 ext

# max V


# min V
