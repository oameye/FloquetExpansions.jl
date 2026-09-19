Reference-only implementation note for issue #169/#98 integration.

The open-leg evaluation plan is a graded extension of the sparse Bloch evaluation schedule from PR #144. Internal Floquet mismatch controls P/Q projection and homological inversion. Symbolic external output paths are interned separately and survive P projection. Vertex order is explicit so A1 and A2 may contribute at different powers of delta. This file is not production API.
