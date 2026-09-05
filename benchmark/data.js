window.BENCHMARK_DATA = {
  "lastUpdate": 1788546264125,
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
          "id": "7704fc2fff94aec6385fba8a3797fd60dbe9b50e",
          "message": "refactor: clean new implementation starting from SQAv0.10 (#37)\n\n* cleanup\n\n* dev docs\n\n* periodic operator\n\n* collector\n\n* engine\n\n* quasi-energy operator\n\n* quasi-energy operator\n\n* add initial documentation for FloquetExpansions.jl\n\n* fix: clarify git policy regarding commits and pushes\n\n* remove dev\n\n* add docs/dev/ to .gitignore\n\n* format\n\n* fix: update Julia version in benchmark workflow to 1\n\n* fix: update action versions in format workflow\n\n* fix docs\n\n* add benchmark\n\n* format\n\n* Jet 0.12",
          "timestamp": "2026-08-28T16:14:37+02:00",
          "tree_id": "9069dda898d6fbe3947a47c5020675b1fd5b7b0d",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/7704fc2fff94aec6385fba8a3797fd60dbe9b50e"
        },
        "date": 1787926884764,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 229151.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=135456\nallocs=2560\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 287521,
            "unit": "ns",
            "extra": "gctime=0\nmemory=257168\nallocs=4280\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 382030,
            "unit": "ns",
            "extra": "gctime=0\nmemory=493136\nallocs=7530\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 435242,
            "unit": "ns",
            "extra": "gctime=0\nmemory=247520\nallocs=4821\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 620655,
            "unit": "ns",
            "extra": "gctime=0\nmemory=593264\nallocs=10087\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 1564989.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2094432\nallocs=33240\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 228104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99200\nallocs=2239\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 412281,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191712\nallocs=4296\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "bc800b3ed115d6b4013f0eb0f18ba0f0982d6fe4",
          "message": "Fix  proper expim collection (#49)\n\n* fix: proper expim collection\n\n* remove PeriodicOperator(QAdd, wd)\n\n* simplify(X::PeriodicOperator)\n\n* more refactor\n\n* _funciton rename\n\n* more refactor\n\n* more a vs b example\n\n* add example\n\n* a vs b example integrated in docs\n\n* fix test env\n\n* fix\n\n* fix docs",
          "timestamp": "2026-09-01T09:55:47+02:00",
          "tree_id": "7bbe876542842aa5d010f2287effa90ac5c8cd45",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/bc800b3ed115d6b4013f0eb0f18ba0f0982d6fe4"
        },
        "date": 1788249819795,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 504863.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=190176\nallocs=3864\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 611002,
            "unit": "ns",
            "extra": "gctime=0\nmemory=316688\nallocs=5592\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 778705,
            "unit": "ns",
            "extra": "gctime=0\nmemory=563248\nallocs=8870\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 992815,
            "unit": "ns",
            "extra": "gctime=0\nmemory=375696\nallocs=7972\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1338525.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=747968\nallocs=13356\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2923127,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2388272\nallocs=36832\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 475418,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152960\nallocs=3533\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 945658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=320592\nallocs=7489\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c060a4770707ef991ca70d25672ae67f54ac289a",
          "message": "feat: Liouvillian van vleck (#50)\n\n* liouvillian\n\n* liouviliian van vleck\n\n* small style changes\n\n* cleanup\n\n* more review changes\n\n* implement review\n\n* change jump API\n\n* fix dev docs\n\n* small architectural refactor\n\n* more restructure\n\n* dissaptive Quasinerhy once more\n\n* remove bad api\n\n* actions(L::Liouvillian)\n\n* better agent docs\n\n* implement docs\n\n* impement docs\n\n* format\n\n* proper manuel\n\n* first iteration theory\n\n* fix docs\n\n* set draft to false in documentation build",
          "timestamp": "2026-09-03T12:50:40+02:00",
          "tree_id": "6072145d6c0b0ddf6e759baf3da5f60d8c253dfd",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/c060a4770707ef991ca70d25672ae67f54ac289a"
        },
        "date": 1788433132855,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 502707,
            "unit": "ns",
            "extra": "gctime=0\nmemory=193664\nallocs=3901\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 607428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330720\nallocs=5713\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 779260.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596416\nallocs=9165\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 963889.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=379184\nallocs=8013\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1335631,
            "unit": "ns",
            "extra": "gctime=0\nmemory=781040\nallocs=13587\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2931080.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2458240\nallocs=37511\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 477295,
            "unit": "ns",
            "extra": "gctime=0\nmemory=154336\nallocs=3554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 921928,
            "unit": "ns",
            "extra": "gctime=0\nmemory=321968\nallocs=7514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4851d72d3e18a77b2bca7d025568f4f5886123b7",
          "message": "chore: remove unused source entries for SecondQuantizedAlgebra in Project.toml files (#53)",
          "timestamp": "2026-09-03T16:25:42+02:00",
          "tree_id": "688624d042a8329f6adac6736fb62fc683230ae4",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/4851d72d3e18a77b2bca7d025568f4f5886123b7"
        },
        "date": 1788446002365,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 336467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=193664\nallocs=3901\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 423419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330720\nallocs=5713\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 581947,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596416\nallocs=9165\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 682161.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=379184\nallocs=8013\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 986704,
            "unit": "ns",
            "extra": "gctime=0\nmemory=778992\nallocs=13555\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2365649,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2434432\nallocs=37143\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 317349.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=154336\nallocs=3554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 653043,
            "unit": "ns",
            "extra": "gctime=0\nmemory=321968\nallocs=7514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a808693d0ac8f454d602163bae5b68ab93d3d615",
          "message": "feat: add exact GKSL coordinate extraction (#66)\n\n* feat: add dissipative frame and GKSL coordinates\n\n* feat: export GKSL coordinate API\n\n* test: cover exact GKSL coordinate extraction\n\n* test: compare Hamiltonian modulo identity in operator gauge\n\n* test: fix isolated GKSL test setup\n\n* fix: specialize Floquet GKSL accessors\n\n* fix: include specialized GKSL Floquet accessors\n\n* docs: add GKSL coordinate API reference\n\n* docs: add GKSL coordinates to manual navigation\n\n* test: use structural zero assertions\n\n* ci: capture JuliaFormatter output\n\n* style: format GKSL coordinate implementation\n\n* chore: remove temporary formatter capture workflow\n\n* docs: keep GKSL coordinates in existing manual pages\n\n* docs: fold GKSL coordinates into existing manual\n\n* docs: integrate dissipative coordinates into system manual\n\n* docs: integrate Floquet GKSL accessors into expansion manual\n\n* ci: preserve successful documentation caches\n\n* fix: align GKSL coordinates with specification\n\n* docs: centralize Floquet GKSL accessors\n\n* api: export GKSL coordinate error\n\n* test: cover coordinate errors and higher-order extraction\n\n* style: format GKSL exports\n\n* style: format GKSL coordinate errors\n\n* style: wrap coordinate validation error\n\n* ci: capture exact formatter diff\n\n* ci: make formatter capture self-contained\n\n* style: apply exact JuliaFormatter output\n\n* ci: remove formatter capture helper\n\n* refine GKSL coordinate errors and frame invariants\n\n* keep GKSL coordinate errors internal\n\n* test GKSL coordinates across operator algebras\n\n* preserve inference for immutable dissipative frames\n\n* avoid coefficient matrix internals in public API tests\n\n* ci: capture formatter output\n\n* format GKSL coordinate implementation\n\n* format GKSL exports\n\n* ci: remove formatter capture helper\n\n* make dissipative frame construction type stable",
          "timestamp": "2026-09-04T20:16:25+02:00",
          "tree_id": "b9e697a0e1cae16180c75ccaf2bab79ae5f9a2a6",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/a808693d0ac8f454d602163bae5b68ab93d3d615"
        },
        "date": 1788546262863,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 487117,
            "unit": "ns",
            "extra": "gctime=0\nmemory=193664\nallocs=3901\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 591210,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330720\nallocs=5713\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 766467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596416\nallocs=9165\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 969470,
            "unit": "ns",
            "extra": "gctime=0\nmemory=379184\nallocs=8013\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1313811,
            "unit": "ns",
            "extra": "gctime=0\nmemory=778992\nallocs=13555\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2870786,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2434432\nallocs=37143\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 469665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=154336\nallocs=3554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 937926,
            "unit": "ns",
            "extra": "gctime=0\nmemory=321968\nallocs=7514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}