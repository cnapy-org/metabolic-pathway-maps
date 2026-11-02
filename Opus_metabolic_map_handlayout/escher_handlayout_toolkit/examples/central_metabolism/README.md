# Example: iML1515 central metabolism (BRENDA draft, 91 reactions)

```bash
escher-handlayout run spec.py --out out --render
```

Reproduces `reference_v7.json`: identical node positions, curves and reaction labels. The canvas
differs by a few tens of units because exit label 1155 is now assigned by the automatic owner rule.

**Layout:**
- **Spine:** glycolysis spine → pyruvate → acetyl-CoA → TCA ring (R = 1050, glyoxylate shunt inside).
- **Right side:** PPP beside upper glycolysis; ED mirrored above 6-phosphogluconate; fermentations
  next to pyruvate and acetyl-CoA.
- **Left side:** gluconeogenesis/anaplerosis.
- **Ring families:** malate↔OAA (5 reactions) and fumarate↔succinate (3) use one concentric lane per
  reaction, with cofactors in the gaps.
- **Merges:** 10 key metabolites. Details are in the spec comments.

**Known weak spots** that the QA report shows for this example and that were accepted in the session:
- The PYK1–4 stacked rows have neighbouring cofactor dots about 30 apart. The lines of the stack run
  through its own rows.
- The ETHAAL nh4 dot sits next to the ALCD2x nadh dot (about 28 apart).
- The LDH_D line passes 11 units from PYK4's ctp dot.
- The L_LACD2/3 stack shows one stack-internal crossing.

Fixing these would be a good first exercise with the QA loop: for example
`CURVES['bend_override']`, PYK row spacing 160, or a `sides=` override for nh4.
