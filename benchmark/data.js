window.BENCHMARK_DATA = {
  "lastUpdate": 1790762985783,
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
          "id": "6e45ee6131429599744e6cb5cbaf61aafee6f145",
          "message": "feat: add completion state, microscopic provenance, and strict jump-rate semantics (#67)\n\n* feat: add completion state, microscopic provenance, and strict jump-rate semantics\n\n* ci: apply targeted #57 CI corrections\n\n* ci: retrigger targeted #57 corrections\n\n* ci: run targeted #57 corrections on branch push\n\n* ci: stage #57 correction script\n\n* ci: simplify targeted #57 correction runner\n\n* fix: address #57 CI regressions\n\n* ci: remove temporary #57 correction workflow\n\n* ci: remove temporary #57 correction script\n\n* ci: apply exact formatter output for #57\n\n* style: format #57 engine changes\n\n* fix: restore micromotion inference\n\n* fix: enforce jump-rate semantics without symbolic complex simplification\n\n* ci: remove temporary formatter workflow\n\n* ci: add targeted #57 repair script\n\n* ci: run targeted #57 repair\n\n* ci: retarget #57 repair to current source\n\n* ci: rerun targeted #57 repair\n\n* ci: narrow #57 repair to remaining gaps\n\n* ci: trigger narrowed #57 repair\n\n* ci: make #57 repair whitespace-robust\n\n* ci: retry robust #57 repair\n\n* fix: close remaining #57 CI gaps\n\n* ci: remove one-shot #57 repair workflow\n\n* ci: remove one-shot #57 repair script\n\n* ci: apply final #57 formatter pass\n\n* style: format final #57 fixes\n\n* ci: remove final #57 formatter helper\n\n* fix: reject complex symbolic jump rates\n\n* docs: separate positive completion manual\n\n* docs: add positive completion manual\n\n* docs: add completion manual to navigation\n\n* docs: keep completion manual focused\n\n* docs: rebalance completion manual and API docstrings\n\n* test: assert public completion guidance\n\n* fix: validate compact symbolic jump phases\n\n* fix completion dispatch against architecture spec\n\n* fix completion order convention in manual\n\n* simplify physical channel lowering and jump validation\n\n* tighten physical Floquet constructor semantics\n\n* expand public completion-state boundary tests\n\n* cover physical jump-rate boundary cases\n\n* ci: temporarily format #57 source changes\n\n* style: format #57 source changes\n\n* refactor completion helper names\n\n* ci: apply #57 internal naming cleanup\n\n* refactor: clean #57 internal helper names\n\n* ci: remove temporary #57 formatter workflow\n\n* ci: apply final #57 test fixes\n\n* fix: preserve real symbolic rates and semantic channel test\n\n* ci: remove temporary #57 fix workflow\n\n* test: compare collapse and rate channels in GKSL coordinates\n\n---------\n\nCo-authored-by: github-actions[bot] <41898282+github-actions[bot]@users.noreply.github.com>",
          "timestamp": "2026-09-05T08:49:23+02:00",
          "tree_id": "e7f1968a88b69f0fa6494ea110314ff1c1cd489f",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/6e45ee6131429599744e6cb5cbaf61aafee6f145"
        },
        "date": 1788591423331,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 336293,
            "unit": "ns",
            "extra": "gctime=0\nmemory=193552\nallocs=3899\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 424491,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330656\nallocs=5712\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 580888.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596400\nallocs=9165\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 678738,
            "unit": "ns",
            "extra": "gctime=0\nmemory=379072\nallocs=8011\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 980276,
            "unit": "ns",
            "extra": "gctime=0\nmemory=780976\nallocs=13586\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2381588,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2458224\nallocs=37511\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 316669.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=154336\nallocs=3554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 650202,
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
          "id": "59d1b5dea4e9d8ace3f873ebd971b971b8ce8c02",
          "message": "docs: add CP-preserving completion theory and user guide (#81)\n\n* docs: add CP-completion theory and user guide\n\n* refactor: tighten completion API and benchmark taxonomy\n\n* docs: align CP completion with manual theory and Literate structure\n\n* refactor: clarify Floquet completion API contracts\n\n* refactor: remove redundant coefficient assertions\n\n* test: import qualified completion expert API\n\n* docs: fix completion references and document order\n\n* docs: separate Floquet theory from completion API\n\n* ci: expose formatter diff temporarily\n\n* ci: apply formatter output temporarily\n\n* style: apply JuliaFormatter\n\n* ci: restore formatter check\n\n* ci: apply focused cleanup temporarily\n\n* ci: tighten focused cleanup checks\n\n* refactor: remove superseded completion helpers\n\n* ci: restore formatter check\n\n* docs: strengthen completion algorithm docstrings\n\n* docs: keep completion manual workflow focused\n\n* docs: keep CP completion theory implementation free\n\n* docs: remove meta commentary from CP example\n\n* docs: separate Floquet theory from package API\n\n* docs: remove package commentary from high-frequency theory\n\n* style: restore formatted completion type source\n\n* docs: clarify spectral completion contract\n\n* docs: sharpen CP completion theory\n\n* docs: label Floquet theory references\n\n* docs: label high-frequency theory\n\n* docs: surface positive completion on landing page\n\n* docs: make active Gram matching perturbative\n\n* docs: qualify automatic completion frames\n\n* docs: fix generic Floquet generator notation\n\n* docs: remove completion implementation caveat\n\n---------\n\nCo-authored-by: github-actions[bot] <41898282+github-actions[bot]@users.noreply.github.com>",
          "timestamp": "2026-09-09T10:08:27+02:00",
          "tree_id": "4d40c5e8992bc00cd8add10851d1655ace5bffd0",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/59d1b5dea4e9d8ace3f873ebd971b971b8ce8c02"
        },
        "date": 1788941876638,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 450828,
            "unit": "ns",
            "extra": "gctime=0\nmemory=193456\nallocs=3896\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 576522.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330528\nallocs=5710\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 749900.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596240\nallocs=9164\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 925839,
            "unit": "ns",
            "extra": "gctime=0\nmemory=378976\nallocs=8008\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1302107.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=778800\nallocs=13552\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2944007,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2426064\nallocs=37014\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 429410,
            "unit": "ns",
            "extra": "gctime=0\nmemory=154336\nallocs=3554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 888734,
            "unit": "ns",
            "extra": "gctime=0\nmemory=321968\nallocs=7514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 11075350,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3663192\nallocs=85107\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 6239814,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2532192\nallocs=50954\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 4029112,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2239344\nallocs=29478\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 3925711,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1930736\nallocs=27281\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 19056068,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5956440\nallocs=149981\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "a0ce2447f736943095e7a8e7bac277c309beab49",
          "message": "ci: add CodeRatchet quality gate (#83)",
          "timestamp": "2026-09-10T20:24:38+02:00",
          "tree_id": "c4b2724e72615c707e7e26d838f3cb936e0f5af5",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/a0ce2447f736943095e7a8e7bac277c309beab49"
        },
        "date": 1789065166853,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 422707.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191920\nallocs=3868\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 532451,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329632\nallocs=5702\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 721121,
            "unit": "ns",
            "extra": "gctime=0\nmemory=597392\nallocs=9220\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 869212,
            "unit": "ns",
            "extra": "gctime=0\nmemory=376384\nallocs=7996\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1241265,
            "unit": "ns",
            "extra": "gctime=0\nmemory=782352\nallocs=13732\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2964006,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2467312\nallocs=38372\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 403139,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152928\nallocs=3530\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 845026,
            "unit": "ns",
            "extra": "gctime=0\nmemory=319504\nallocs=7506\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 8761878,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2953056\nallocs=74712\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 5455797,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2057528\nallocs=48238\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 2765822,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1127944\nallocs=18879\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2668257.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=939536\nallocs=17591\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 14334957,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4740096\nallocs=124559\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "24815fdca283ec12f7016b6dfc5c6dc96e18c0bd",
          "message": "docs: show methods (#87)\n\n* docs: show methods\n\n* format",
          "timestamp": "2026-09-10T20:28:42+02:00",
          "tree_id": "2a9dd2d657234822ab08813cb32feb5e6edfe67b",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/24815fdca283ec12f7016b6dfc5c6dc96e18c0bd"
        },
        "date": 1789065217363,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 520407.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191920\nallocs=3868\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 629725.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329632\nallocs=5702\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 803694,
            "unit": "ns",
            "extra": "gctime=0\nmemory=597392\nallocs=9220\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 1015527,
            "unit": "ns",
            "extra": "gctime=0\nmemory=376384\nallocs=7996\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1368599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=782352\nallocs=13732\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2975299,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2467312\nallocs=38372\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 500930,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152928\nallocs=3530\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 998456,
            "unit": "ns",
            "extra": "gctime=0\nmemory=319504\nallocs=7506\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 9853682,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2953056\nallocs=74712\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 6199066.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2057528\nallocs=48238\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 3096455,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1127944\nallocs=18879\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 3009394,
            "unit": "ns",
            "extra": "gctime=0\nmemory=939536\nallocs=17591\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 15888405.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4740096\nallocs=124559\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "a5cd0170a54367bff0994980368e1aca2de7a258",
          "message": "docs: rework the completion docs and add rich Liouvillian displays (#116)\n\n* better docs\n\n* add Standars.md back\n\n* fix: remove machine-local paths from the test and docs environments\n\n`test/Project.toml` pinned FloquetExpansions to an absolute path under a developer home\ndirectory, so every CI test job failed at instantiation with \"expected package\nFloquetExpansions to exist at path ...\". The entry is redundant in any case: Pkg.test\ndevelops the package under test on its own.\n\n`docs/Project.toml` pointed SecondQuantizedAlgebra at a local clone, so the docs build\nfailed with \"Path /home/... does not exist\". It now uses the same GitHub URL as the root\nproject.\n\nPkg writes both entries when it resolves inside a local checkout. Neither belongs in a\ncommitted environment.\n\nAlso corrects a misspelling that the typos gate rejected.\n\n* chore: record the argument count of the rich display methods\n\nAdding `show(io, ::MIME\"text/plain\", x)` to the GKSL-coordinate and Gram-completion files\nraised their maximum argument count from two to three, which the complexity ratchet held.\nThe count cannot be lowered, because the three-argument form is Base's signature for a\nrich display method.\n\nThe baseline is measured rather than authored, so it is refreshed wholesale. No capped\nmetric moved: every `arg_over`, `cyc_over` and `cog_over` stays at zero, and only the\nper-file sums and those two maxima changed.\n\n* fix: build channel LaTeX fragments through the display path\n\n`latex_fragment` called `Latexify.latexify(x; env=:raw)`. Passing an environment routes\nthrough `_latexraw`, which has no method for the symbolic type backing a `Coeff` on the\nversions that resolve for Julia 1.10, so every channel LaTeX display threw\n\"cannot latexify objects of type ...\". Julia 1.13 resolved versions where it happened to\nwork, which is why local runs stayed green.\n\nRender through `show(::MIME\"text/latex\")`, the path both SecondQuantizedAlgebra and\nSymbolics support, and strip the outermost delimiter pair from the result. The stripped\noutput is identical on 1.10 and 1.13, and Latexify is no longer a dependency.\n\n* fix: bound SymbolicUtils below the 4.46.5 coefficient regression\n\nSymbolicUtils 4.46.5 changed `poly_to_gcd_form` so that a polynomial coefficient set mixing\nan exact `Rational` with a complex whose imaginary part vanishes takes the rational branch:\n`isinteger(4.0 + 0.0im)` is true, so `all_rat` stays set even though `any_complex` is also\nset. That branch calls `rationalize`, which returns `Complex{Rational{Int64}}`, and then\n`numerator` on it.\n\n`Base.numerator(::Complex{<:Rational})` exists only from Julia 1.13, so the regression\nsurfaces as a MethodError on 1.10 and 1.11 and is merely hidden on 1.13. Spectral completion\nproduces exactly that coefficient mix, since the retained rates are exact rationals while the\nbranch square roots are floats.\n\nBisected: 4.46.4 passes, 4.46.5 fails. SymbolicUtils becomes a direct dependency because Pkg\nignores a compat entry for a package that is not one. Lift the bound once this is fixed\nupstream.\n\n* fix: bound SymbolicUtils in the test environment, not the package\n\nReworks the previous commit. A package-level bound is unsatisfiable: the documentation\nenvironment pins QuantumCumulants to master, which requires SymbolicUtils 4.46.5 or newer,\nso capping the package at 4.46.4 leaves the docs build with no resolvable version.\n\nThe bound therefore lives in the test environment, which is where the failing gate runs, and\nSymbolicUtils is no longer a direct dependency of the package. The docs build resolves freely\nand is unaffected in practice, because it runs on Julia 1.13 where the regression is masked\nby the Base method added there.\n\nThe underlying defect is upstream and still wants reporting: `poly_to_gcd_form` treats a\ncoefficient set mixing an exact Rational with a complex whose imaginary part vanishes as\nrational, because `isinteger(4.0 + 0.0im)` is true.",
          "timestamp": "2026-09-14T10:18:19+02:00",
          "tree_id": "723c100e8fee536236ab91fba72e6fc9da231949",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/a5cd0170a54367bff0994980368e1aca2de7a258"
        },
        "date": 1789374333078,
        "tool": "julia",
        "benches": [
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 513118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=190000\nallocs=3828\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 613882.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=327616\nallocs=5660\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 786342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=595280\nallocs=9176\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 995134,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372256\nallocs=7910\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1350455,
            "unit": "ns",
            "extra": "gctime=0\nmemory=780432\nallocs=13684\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 3037985,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2479120\nallocs=38562\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 491082,
            "unit": "ns",
            "extra": "gctime=0\nmemory=151008\nallocs=3490\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 970407.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=315376\nallocs=7420\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 9678620,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2943248\nallocs=74497\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 5974822,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2048584\nallocs=48037\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 3050393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1124728\nallocs=18812\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2964514.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=936320\nallocs=17524\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 15530251,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4729056\nallocs=124317\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "c83df19e817193a57fd53eb27a8379d1a5a8ea69",
          "message": "refactor: land generic HD/BF Floquet HFE core (#323)\n\n* refactor: land generic Hori-Deprit and Bloch-Feshbach HFE core\n\n* Add documentation for expansion algorithms and Bloch/Feshbach projection\n\n- Introduced a new manual page for expansion algorithms, detailing the selection process and available algorithms.\n- Updated the Floquet expansion documentation to clarify the separation of gauge and expansion algorithm choices.\n- Added references for Bloch/Feshbach projection theory, including detailed explanations of the wave operator and its relation to effective generators.\n- Enhanced the `ExpansionAlgorithm` type with descriptions for `HoriDeprit` and `BlochFeshbach`, including usage examples.\n- Improved comments in the engine and periodic operator files to clarify the role of gauge and algorithm in Floquet expansions.\n\n* restructure\n\n* Refactor Bloch-Feshbach implementation and introduce Van Vleck words\n\n- Removed unused files related to connected log and van Vleck normalization.\n- Introduced new structures and functions for handling Lyndon words and brackets.\n- Updated the Floquet expansion implementation to utilize the new Lie series plan.\n- Enhanced the wave operator logic to support the new Van Vleck words structure.\n- Improved the overall organization of the Bloch-Feshbach module for better clarity and maintainability.\n- Adjusted the completion and gauge handling to reflect changes in the underlying structures.\n\n* refactor: update file paths and provenance in baseline TOML files; enhance documentation and tests for PeriodicGenerator\n\n* fix spelling\n\n* final review\n\n* refactor: organize and update expansion algorithm exports in tests\n\n* fix",
          "timestamp": "2026-09-26T23:21:50+02:00",
          "tree_id": "2eacfbf87e7a5e82ee89ecb58958419220792b2b",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/c83df19e817193a57fd53eb27a8379d1a5a8ea69"
        },
        "date": 1790458301517,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1487434,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1297216\nallocs=13896\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 14854643.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18950896\nallocs=180281\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 5012500.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5975048\nallocs=58078\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 277437584,
            "unit": "ns",
            "extra": "gctime=38574929\nmemory=225354016\nallocs=2551696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 591506,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314176\nallocs=5166\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1108602,
            "unit": "ns",
            "extra": "gctime=0\nmemory=983448\nallocs=12253\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 728478.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=554336\nallocs=8560\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 2035950,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2380000\nallocs=33972\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1306635,
            "unit": "ns",
            "extra": "gctime=0\nmemory=780256\nallocs=12496\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 7996043.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8074600\nallocs=90671\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2837719,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2398512\nallocs=37336\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 51727906,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43946368\nallocs=578681\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 522311,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191024\nallocs=3852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 623892.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328640\nallocs=5684\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 791898,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596304\nallocs=9200\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 1020913,
            "unit": "ns",
            "extra": "gctime=0\nmemory=374368\nallocs=7962\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1357209,
            "unit": "ns",
            "extra": "gctime=0\nmemory=780240\nallocs=13696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2913152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2465104\nallocs=38334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 497134,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152032\nallocs=3514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 991923,
            "unit": "ns",
            "extra": "gctime=0\nmemory=317488\nallocs=7472\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 7948064.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2348992\nallocs=58102\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 5335237,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1797928\nallocs=40936\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 2429974,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1060008\nallocs=16903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2351647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=871600\nallocs=15615\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 13041511,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3897152\nallocs=101334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "8fe4108704b574e347a13cdfc7fb21eb300c9b64",
          "message": "fix: complete coherent Liouvillians in an empty dissipative frame (#335)\n\npositive_completion(expansion, Gram()) and positive_completion(expansion,\nSpectral()) threw an ArgumentError from automatic frame discovery when the\nexpansion had no dissipative directions, for example a purely coherent\nLiouvillian. A coherent generator is trivially GKSL: its completion has an\nempty dissipative frame, no channels, a 0x0 Kossakowski matrix, and the\nHamiltonian action of the retained coherent Hamiltonian as its effective\ngenerator.\n\nAutomatic frame discovery now assembles its frame without the nonempty guard\nof the public DissipativeFrame constructor. DissipativeFrame() still rejects\nan empty user frame, and GKSL extraction in the empty frame still rejects any\ndissipative residual, so a dissipative expansion cannot reach this path.",
          "timestamp": "2026-09-27T07:38:26+02:00",
          "tree_id": "60dce869880990f8bf53b83e8fe6b1e866568a96",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/8fe4108704b574e347a13cdfc7fb21eb300c9b64"
        },
        "date": 1790487851469,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1570651,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1297216\nallocs=13896\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 15634636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18963904\nallocs=180245\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 5396992,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5976072\nallocs=58102\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 296292884,
            "unit": "ns",
            "extra": "gctime=42180517\nmemory=223825232\nallocs=2550560\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 617163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314176\nallocs=5166\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1163302,
            "unit": "ns",
            "extra": "gctime=0\nmemory=983448\nallocs=12253\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 761599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=554336\nallocs=8560\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 2142959,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2380000\nallocs=33972\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1358710,
            "unit": "ns",
            "extra": "gctime=0\nmemory=780256\nallocs=12496\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 8354765,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8074600\nallocs=90671\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2916201,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2398512\nallocs=37336\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 47175710.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43946368\nallocs=578681\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 538145,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191024\nallocs=3852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 631801,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328640\nallocs=5684\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 809390.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596304\nallocs=9200\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 1070287.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=374368\nallocs=7962\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1409927,
            "unit": "ns",
            "extra": "gctime=0\nmemory=780240\nallocs=13696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2992821,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2465104\nallocs=38334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 508338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152032\nallocs=3514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 1008696.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=317488\nallocs=7472\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 8155252.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2348992\nallocs=58102\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 5552845,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1797928\nallocs=40936\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 2492719,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1060008\nallocs=16903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2402056,
            "unit": "ns",
            "extra": "gctime=0\nmemory=871600\nallocs=15615\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 13349099,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3897152\nallocs=101334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "27313f939663857083c1e6aa18f2066c2b540f09",
          "message": "refactor: restructure src and test into layer folders (#344)\n\n* refactor: restructure src and test into layer folders\n\nSplit the flat source tree into generators/, words/, expansion/, gksl/ and\ncompletion/, included in that order. Tests mirror the same folders.\n\nThree types sat in a later layer than the code that needed them. The\nprovenance types move beside the channels that build them, Completion and\nUncompleted move into expansion/, and antiderivative(G, ::VanVleck) moves\ninto gauges.jl. The accessors shared by Gram and Spectral move out of the\nGram file into completion/accessors.jl.\n\nliouvillian.jl splits into the map algebra and channels.jl.\ngksl_coordinates.jl splits into dissipative_frame.jl and coordinates.jl.\nNo identifier, signature or behaviour changes: every named binding and\nmethod signature matches the pre-move module.\n\narchitecture.md gains a layering rule and a per-file ownership table.\nADR test pointers follow the moved tests. Ratchet baselines are refreshed;\ntheir totals equal a fresh measurement of main.\n\n* docs: drop the restructure plan",
          "timestamp": "2026-09-27T14:42:11+02:00",
          "tree_id": "23d5e007869a68d5776e953538de087e86216593",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/27313f939663857083c1e6aa18f2066c2b540f09"
        },
        "date": 1790513293849,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1136671,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1297216\nallocs=13896\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 13346243,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18946056\nallocs=180253\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 4346516.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5995800\nallocs=58056\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 254568744,
            "unit": "ns",
            "extra": "gctime=41427443\nmemory=223024752\nallocs=2551042\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 409299,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314176\nallocs=5166\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 859058.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=983448\nallocs=12253\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 536287.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=554336\nallocs=8560\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1650586,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2380000\nallocs=33972\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 948771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=780256\nallocs=12496\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 6703478,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8070440\nallocs=90661\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2283698,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2434480\nallocs=37956\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 39726178,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43917632\nallocs=576647\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 343931.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191024\nallocs=3852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 433670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328640\nallocs=5684\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 583297,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596304\nallocs=9200\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 705900,
            "unit": "ns",
            "extra": "gctime=0\nmemory=374368\nallocs=7962\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 992045,
            "unit": "ns",
            "extra": "gctime=0\nmemory=782544\nallocs=13736\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2346476,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2501072\nallocs=38954\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 325584.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152032\nallocs=3514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 677237,
            "unit": "ns",
            "extra": "gctime=0\nmemory=317488\nallocs=7472\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 5628824,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2348992\nallocs=58102\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 3773901,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1797928\nallocs=40936\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 2427041,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1220328\nallocs=21241\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2335835,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1031920\nallocs=19953\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 9398498.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3897152\nallocs=101334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "195b75effac470a4f726233c2599c237be380d2a",
          "message": "fix: seed automatic frames with static channel harmonics (#349)\n\n* fix: seed automatic frames with static channel harmonics\n\nAutomatic dissipative-frame discovery appended the microscopic channel\noperators as frame directions without Fourier lowering them, so a\ntime-dependent jump or collapse operator entered the frame with its\ndrive phases and Gram completion failed on a non-Hermitian retained\nKossakowski series.\n\nThe provenance now records, in user channel order, the Fourier harmonics\nof each channel operator with the time average first, and discovery\nseeds the frame from those. A static channel still seeds its own\noperator. The seed-reference enum and struct are removed, since\ndiscovery was their only consumer.\n\n* Remove CHANGELOG.md\n\n* refactor: pass lowered harmonics to append_frame_seeds!\n\nThe Ratchet complexity gate holds src/generators/channels.jl at a\nmaximum of three arguments per function. append_frame_seeds! took the\noperator, drive frequency and time separately; it now takes the lowered\nPeriodicGenerator, and microscopic_provenance performs the lowering at\nthe call site. Behaviour is unchanged.",
          "timestamp": "2026-09-28T10:52:16+02:00",
          "tree_id": "a93232868fc43b113c0e4b035445f84741d43527",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/195b75effac470a4f726233c2599c237be380d2a"
        },
        "date": 1790585919502,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1475127,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1311200\nallocs=14090\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 16673249,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18949088\nallocs=180444\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 5510088,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5989032\nallocs=58272\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 316118567,
            "unit": "ns",
            "extra": "gctime=45979457\nmemory=225368128\nallocs=2551958\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 520151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314176\nallocs=5166\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1086557.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=983448\nallocs=12253\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 689390.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=554336\nallocs=8560\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 2110721,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2380000\nallocs=33972\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1207525,
            "unit": "ns",
            "extra": "gctime=0\nmemory=779680\nallocs=12486\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 8688157,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8066568\nallocs=90596\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2925591,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2405424\nallocs=37456\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 55843764,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44047040\nallocs=579801\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 443667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191024\nallocs=3852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 555690,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328640\nallocs=5684\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 746081,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596304\nallocs=9200\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 913172.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=374368\nallocs=7962\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1270464,
            "unit": "ns",
            "extra": "gctime=0\nmemory=782544\nallocs=13736\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2999098,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2472016\nallocs=38454\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 419360,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152032\nallocs=3514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 878140,
            "unit": "ns",
            "extra": "gctime=0\nmemory=317488\nallocs=7472\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 7284637,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2348992\nallocs=58102\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 4910882.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1797928\nallocs=40936\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 3119034,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1220328\nallocs=21241\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2989833.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1031920\nallocs=19953\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 12263772,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3897152\nallocs=101334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "b8f00a87897d063b8049f04a2ab62c61a10f7b82",
          "message": "fix: canonical Liouvillian keys; docs: recurrence-cost plot (#345)\n\n* docs: plot the Lie-transform and Bloch recurrence operation counts\n\nCompare the N(N^2-1)/3 commutators of the Hori-Deprit triangle with the\n(N-1)(N+2)/2 products of the Bloch recurrence, and state that the counts\nexclude harmonic convolution and the van Vleck normalization.\n\n* fix: key Liouvillian terms by canonical monomial pairs\n\nTerms were keyed by raw (A, B) factor pairs, so equal maps with different\nfactorizations never merged and a zero map could hold terms. The\ncommutator of hamiltonian_action(a'a) and dissipator(a) kept four terms,\nand a test asserted it was nonzero.\n\nadd_term! now splits each factor into unit monomials after the NLevel\ncompleteness relation, so iszero and == are map equality. Factors with\nbound symbolic sums keep the whole-factor key. On a driven dissipative\nqubit at order 6 this cuts effective terms from 161 to 8 and Hori-Deprit\ntime from 2.6 s to 0.63 s. ADR 0007 is amended.\n\n* refactor: split canonical Liouvillian term insertion into helpers\n\nKeeps add_term! within the complexity ratchet.\n\n* add scaling plot\n\n* add sclaing plot",
          "timestamp": "2026-09-30T11:59:52+02:00",
          "tree_id": "0f27014dc15d0798209cefed11a6dab16a858b5c",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/b8f00a87897d063b8049f04a2ab62c61a10f7b82"
        },
        "date": 1790762982483,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1785890,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2381264\nallocs=21223\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 6195653,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12657624\nallocs=94780\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 4457457,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9266672\nallocs=70082\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 54352869,
            "unit": "ns",
            "extra": "gctime=7737497\nmemory=88203104\nallocs=635490\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 606193,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314176\nallocs=5166\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1159478,
            "unit": "ns",
            "extra": "gctime=0\nmemory=983448\nallocs=12253\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 752077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=554336\nallocs=8560\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 2082840,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2380000\nallocs=33972\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1329266,
            "unit": "ns",
            "extra": "gctime=0\nmemory=779104\nallocs=12476\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 8146438.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8064872\nallocs=90521\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2808386,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2359088\nallocs=36656\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 50636772,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43571008\nallocs=572687\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 533458,
            "unit": "ns",
            "extra": "gctime=0\nmemory=191024\nallocs=3852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 634161,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328640\nallocs=5684\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 809714,
            "unit": "ns",
            "extra": "gctime=0\nmemory=596304\nallocs=9200\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 1037771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=374368\nallocs=7962\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1401771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=780240\nallocs=13696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2924619,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2425680\nallocs=37654\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 509137,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152032\nallocs=3514\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 1015058,
            "unit": "ns",
            "extra": "gctime=0\nmemory=317488\nallocs=7472\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 8217859,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2579680\nallocs=60038\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 5424962,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2024248\nallocs=42327\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 3808246,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1962864\nallocs=26036\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 3694133,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1670904\nallocs=24260\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 13323416,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128864\nallocs=103263\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}