window.BENCHMARK_DATA = {
  "lastUpdate": 1753115597008,
  "repoUrl": "https://github.com/oameye/FloquetExpansions.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2ff88ea29270c60d6bd1e6bc7bda613bb129cfcb",
          "message": "refactor: fourier components (#5)",
          "timestamp": "2025-07-21T18:09:41+02:00",
          "tree_id": "d662d6c845c3ff1bb7e3bf9555f6167f116a5444",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/2ff88ea29270c60d6bd1e6bc7bda613bb129cfcb"
        },
        "date": 1753114428366,
        "tool": "julia",
        "benches": [
          {
            "name": "Construction/Quasienergy operator/Kerr resonator",
            "value": 103785,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59616\nallocs=1171\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Operation/Rotate/Kerr resonator",
            "value": 136605792.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28243408\nallocs=690021\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Floquet Expansion/Rotating Wave Approximation/Kerr resonator",
            "value": 136342818,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28303008\nallocs=691192\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "distinct": true,
          "id": "c11a566c9035cbfd7571e65ab87c62f5227d677b",
          "message": "docs: add kerr resonator example",
          "timestamp": "2025-07-21T18:31:03+02:00",
          "tree_id": "cad520f3dd818dac0dc44da37e44f5fa72330c8b",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/c11a566c9035cbfd7571e65ab87c62f5227d677b"
        },
        "date": 1753115595983,
        "tool": "julia",
        "benches": [
          {
            "name": "Construction/Quasienergy operator/Kerr resonator",
            "value": 101658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59616\nallocs=1171\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Operation/Rotate/Kerr resonator",
            "value": 131188390,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28243392\nallocs=690021\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Floquet Expansion/Rotating Wave Approximation/Kerr resonator",
            "value": 131204444,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28303008\nallocs=691192\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}