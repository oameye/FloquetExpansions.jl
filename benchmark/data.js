window.BENCHMARK_DATA = {
  "lastUpdate": 1792153910465,
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
          "id": "f79a884eb0bb4d69954b56e6876c4af351ee69e9",
          "message": "build: require SecondQuantizedAlgebra 0.13 (#393)\n\nSQA main became 0.13.0 on 2026-10-04, outside the 0.11 - 0.12 compat, so a\nfresh resolve of the [sources] entry failed. SQA 0.13 keeps coefficients exact\nby type, so the positive_completion doctest now prints exact rationals.",
          "timestamp": "2026-10-11T15:54:00+02:00",
          "tree_id": "40c6772997150ec5f0fd9510d13fd58ed0358c48",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/f79a884eb0bb4d69954b56e6876c4af351ee69e9"
        },
        "date": 1791727499146,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1996783,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2541040\nallocs=21406\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 7046902,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13573752\nallocs=95773\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 5028688,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9945504\nallocs=70080\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 64178975,
            "unit": "ns",
            "extra": "gctime=11526495\nmemory=94586656\nallocs=635786\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 738398,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329272\nallocs=5435\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1317989,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1033792\nallocs=13268\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 961824,
            "unit": "ns",
            "extra": "gctime=0\nmemory=584936\nallocs=9220\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 2368192,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2507632\nallocs=36838\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1549690,
            "unit": "ns",
            "extra": "gctime=0\nmemory=815792\nallocs=13224\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 9354977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8356336\nallocs=97556\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 3179680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2437624\nallocs=38124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 47914885,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43814744\nallocs=581925\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 652577,
            "unit": "ns",
            "extra": "gctime=0\nmemory=199440\nallocs=4000\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 753024,
            "unit": "ns",
            "extra": "gctime=0\nmemory=344064\nallocs=5977\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 943059,
            "unit": "ns",
            "extra": "gctime=0\nmemory=628184\nallocs=9860\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 1336693,
            "unit": "ns",
            "extra": "gctime=0\nmemory=390208\nallocs=8264\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1619410,
            "unit": "ns",
            "extra": "gctime=0\nmemory=812048\nallocs=14241\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 3277853,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2506072\nallocs=39122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 626208,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158528\nallocs=3662\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 1294373,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330448\nallocs=7774\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 10799265,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2873888\nallocs=66823\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 7100652,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2214808\nallocs=46329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 2372264,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1776416\nallocs=18105\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2257751,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1468992\nallocs=16500\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 17002535,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4533616\nallocs=113000\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "c4ab31753a288ba746f11732412fe1db4604f5b8",
          "message": "feat: GKSLNormalForm, a completely positive native Floquet gauge (#390)\n\n* test: add native GKLS BF orders 0-2 prototype\n\n* test: format native GKLS BF prototype\n\n* ci: expose formatter output for native BF prototype\n\n* ci: capture JuliaFormatter output for prototype\n\n* test: apply JuliaFormatter to native BF prototype\n\n* ci: restore standard format check\n\n* test: add independent native GKLS HD orders 0-2 oracle\n\n* ci: expose HD prototype formatter output\n\n* test: format native GKLS HD prototype\n\n* ci: restore standard format check\n\n* ci: inspect remaining HD formatter diff\n\n* test: add trailing newline to native HD oracle\n\n* ci: restore standard format workflow\n\n* test: add arbitrary-order native GKLS BF recurrence\n\n* ci: print native BF formatter diff\n\n* style: format arbitrary native GKLS BF oracle\n\n* ci: restore format workflow\n\n* ci: run tests on arbitrary native BF branch\n\n* ci: restore tests workflow\n\n* test: certify graded channel activation after birth\n\n* test: add arbitrary-order native GKLS HD recurrence\n\n* test: fix arbitrary HD state comparison\n\n* ci: print arbitrary HD formatter diff\n\n* style: format arbitrary native GKLS HD oracle\n\n* ci: restore format workflow\n\n* refactor: extract shared Gram/Lyapunov kernel\n\n* research: solve native GKLS static affine PSD section\n\n* test: certify finite-Fock native GKLS Kerr\n\n* test: certify localized Kerr static inverse\n\n* test: certify localized Kerr static ladder\n\nSolve the charged Kerr static equation y_t = d_t s_t + J s_{t+1}, including\nloss recycling, with exact Complex{Rational{BigInt}} arithmetic in the first\nnon-bright modulated-loss sector (q_L, q_R) = (-3, -2). Verify [Lbar, S] = Y on\nthe full finite-Fock superoperator, check that the coefficients on every\noverlapping level and physical support are identical across cutoffs 6, 8, 10,\nand cross-check the ladder against an independent exact elimination over the\nwhole Q = -1 sector. That sector is dissipatively gapped at every level, so the\nsolution is unique.\n\n* test: resolvent section for charged native GKLS static slots\n\nThread a homological-inverse seam through the native static step. The default\nkeeps the dense minimum-norm affine PSD section unchanged. A charged section\ncancels the charged dark-dark target exactly via [Lbar, S] = -Pi_{Q!=0} D[pi\nDelta pi], using either a dense LU of ad_Lbar on Q != 0 or the localized Kerr\nrecycling ladder, and leaves the neutral block to the affine PSD solve over the\nneutral gauge algebra.\n\nOn finite-Fock driven Kerr the ladder equals the dense inverse in the full\nrecurrence at cutoffs 3 and 4, BF equals HD through order 4, and the section\ndiffers from minimum-norm only by a static gauge in ker Phi_4 with spectra\nagreeing to O(eps^5). The Hamiltonian branch is unchanged and physical\nresonances are rejected.\n\n* feat: move native GKLS graded-channel state into src/native\n\nAdd the internal native/ layer for the generator-native GKLS expansion, after\ngksl/ and independent of completion/. GradedChannel{T} holds a jump amplitude\nwith integer onset order (grade onset/2), generic over the scalar type, with\nactive-flag extraction, known products K_n^<, Gram coefficients, and\nprefix-complete storage of corrections and births.\n\nThe research BF/HD recurrences now use this state instead of their private\ncopy. Exact Complex{Rational} tests check the Cauchy-product structure, and a\nJET optimizer workload covers the new functions.\n\n* feat: move the native GKLS static slot core into src/native\n\nMove the floating-point Gram geometry from completion/backend to\ngksl/gram_kernel.jl, so native/ no longer reaches into completion/. The\nFeshbach residual, the only part completion uses, moves beside its caller in\ncompletion/backend/linear_algebra.jl.\n\nAdd native_static_solve: one native static slot in Kossakowski coordinates of\na fixed frame (dark quotient, affine PSD gauge section, newborn channels,\nbright tangent lift, reconstruction check), returning a concretely typed\nNativeStaticSolution. The superoperator representation stays with the caller.\nThe research BF/HD recurrences now run their static slot through it.\n\n* feat: generic native static step and dense Liouville representation\n\nAdd the NativeRepresentation interface and the HomologicalInverse seam to\nsrc/native/static_step.jl. An inverse supplies a static gauge that cancels the\ndark target exactly on the sectors where ad_Lbar is regular, plus the gauge\ndirections of the remaining sector, which keeps the affine PSD section.\nNoHomologicalInverse has no regular sector.\n\nAdd DenseLiouvilleRepresentation, the finite-dimensional backend, with its\nstatic gauge algebra built in the constructor rather than in a global cache.\nThe research recurrences and the charged Kerr section now run through these.\n\nFollow the style guide for internals: positional arguments only, no comments\nin src.\n\n* feat: generic native BF/HD recurrence driver\n\nAdd native_recurrence, the arbitrary-order generator-native GKLS recurrence\nfor the BlochFeshbach and HoriDeprit algorithms, over harmonic series of\nrepresentation elements: intrinsic total-kick offset, one native static step\nper order, and graded-channel bookkeeping.\n\nThe research BF and HD recurrences are now thin wrappers around it. The\nindependent HD Lie-transform series and equation-defect checks stay in the\ntests, so the package HD state is still checked against a separate\nimplementation.\n\n* docs: ADR 0011 for the generator-native GKLS expansion\n\n* feat: exact native static slot for the canonical branch\n\nAdd native_exact_static_solve, the exact counterpart of native_static_solve,\ngeneric over the scalar field. It takes the operator and gauge metrics\nexplicitly, since an exact frame is not orthonormal, and returns the weighted\nleast-squares gauge section with a minimum-norm tie-break, rate-weighted\nnewborn channels from an LDL factorization without square roots, and the exact\nno-bright-rotation tangent lift. An indefinite canonical residual is rejected,\nbecause the affine PSD fallback has no exact analogue.\n\nExact Gauss-Jordan elimination replaces the generic backslash, whose abs-based\npivoting is type-unstable on exact rationals. Tests check agreement with the\nfloat core in orthonormal coordinates, covariance under a non-orthonormal\nframe and under a gauge reparametrization, and rejection of an indefinite\nresidual.\n\n* feat: exact, charge-graded and cutoff-free native GKLS representations\n\n- GradedChannel carries a rate weight; the BF/HD recurrence is generic over\n  float and exact scalars (exact series constants, representation-defined\n  coefficient checks).\n- ExactLiouvilleRepresentation: exact dense backend with an explicit Hilbert\n  metric, so Fock and spin operators are exact in a rescaled basis.\n- ChargeGradedInverse: detects the U(1) charge lattice of the averaged\n  generator and the leading channels, applies the resolvent section on\n  regular charged blocks, and keeps the affine section on the rest.\n- OperatorAlgebra and AlgebraicLiouvilleRepresentation: cutoff-free exact\n  algebra of boson, phase-space and finite-level sites with a polynomial\n  charged homological solve.\n- Sparse exact form and gauge metrics.\n\nCutoff-free driven Kerr has no newborn channel and no static gauge through\norder 3, while finite-Fock truncations birth cutoff-localized channels from\norder 2.\n\n* fix: concrete default Hilbert metric for the exact representation\n\n* feat: exact facial-reduction PSD lift and sparse exact slot\n\n- The exact static slot imposes the linear constraints forced by positivity\n  (rows over the kernel of every gauge-fixed principal block) before the\n  canonical cancellation, so a PSD lift is found exactly where least squares\n  is indefinite; an indefinite fixed block proves that no lift exists.\n- Dark projection by low-rank updates, an identity-metric shortcut and sparse\n  trace products replace the dense n^2 x n^2 form metric.\n\nCutoff-free driven Kerr at order 4 now births one exact channel with rate\n72/625, identical for Bloch-Feshbach and Hori-Deprit.\n\n* perf: cached gauge family, closed-form bright lift, sparse facial products\n\n- The recurrence computes the static gauge images and metric once.\n- The bright lift splits off its dark component in closed form and solves\n  only a k x k system for the active block.\n- Facial-reduction products run over nonzero entries.\n- An exact row-cancellation fallback applies only when the canonical\n  least-squares point is indefinite, and indefinite gauge-fixed blocks are\n  reported as a proof that the gauge family has no PSD lift.\n\n* feat: SQA lowering and exact native expansion driver\n\n- SQALowering maps Fock, phase-space, NLevel and Pauli operators of a periodic\n  Liouvillian into the cutoff-free exact algebra, with exact coefficients and\n  numeric parameter values, and lifts results back to SQA expressions.\n- native_expansion_data runs the native BF/HD recurrence with leading channels\n  from the exact rate-weighted factor of the averaged Kossakowski matrix,\n  enlarges the frame on demand and promotes Int128 to BigInt on overflow.\n- Per-order static gauge-degree escalation; NativePositivityError and\n  NativeFrameError separate positivity and frame failures from bugs.\n\nThe native coefficients equal the symbolic Van Vleck coefficients exactly\nwherever the static gauge vanishes (driven Kerr, orders 0-3, BF and HD).\n\n* feat: GKSLNormalForm gauge for completely positive Floquet expansions\n\nfloquet_expansion(system, GKSLNormalForm(; algorithm), order) runs the exact\nnative BF/HD recurrence on the cutoff-free operator algebra and returns a\nFloquetExpansion carrying a NativeRealization: effective components and\nmicromotion of the GKSL normal form, the finite GKSL generator, its\nHamiltonian, rate-weighted channels and Kossakowski data. The gauge equals\nVanVleck wherever every static gauge vanishes, and returns the van Vleck\ncoefficients for Hamiltonian generators.\n\n- Completed-state accessors accept the native realization; positive\n  completion refuses an already GKSL result.\n- Facial reduction runs one function per step, keeping the exact slot\n  inferable.\n- ADR 0011, CONTEXT.md, the manual and the architecture map describe the\n  gauge, the exact positive section and the cutoff-free backend.\n\n* chore: spell GKSL consistently and label truncated-model Kerr tests\n\nRename the remaining GKLS spellings to GKSL, including the ADR 0011 and research note\nfile names. Label the finite-Fock Kerr fixtures as truncated-model tests that carry the\ncutoff artifact, and fix a typo flagged by the spell check.\n\n* fix: compare native results with the floating Van Vleck reference by tolerance\n\nSince the canonical Liouvillian keys, liouvillian(H) folds the native -im factor into raw\ntrigonometric coefficients, and SQA embeds it as ComplexF64. The Van Vleck reference for\nnumeric models is therefore computed in floating point, while the native expansion stays\nexact. Compare the two with a tolerance.\n\nAlso clear the JETLS diagnostics in the native files: unused arguments and locals, an\nunsorted export list, and Bool assertions where inference saw Missing.\n\n* chore: clear the remaining JETLS warnings and refresh the lsp baseline\n\n* fix: keep exact drive amplitudes exact through Liouvillian construction\n\nCanonical term insertion multiplies the -im of a Hamiltonian action into raw trigonometric\ncoefficients, and SQA embeds a native factor of a raw expression as ComplexF64. The 2//5 of\n(2//5) cos(ωt) became 0.4, harmonics returned float amplitudes, and both the Van Vleck\nrecurrence and the native lowering lost exactness: a cos/sin drive gave a garbage newborn\nrate. Promote Gaussian-integer native factors to exact constants before they meet a raw\ncoefficient in term insertion, scaling, and composition.\n\nRestore exact native-versus-Van Vleck comparisons in the native tests and drop a local\npath source that slipped into test/Project.toml.\n\n* feat: spin sites, symbolic parameters, and exact finite assembly for GKSLNormalForm\n\n- Spin operators lower to a size-independent U(su(2)) site in the PBW basis\n  S+^a Sz^b S-^c. Monomial slots are now per site, three for spin and two otherwise.\n- With symbolic parameters, the structure is decided at two generic rational points. If\n  no static gauge is needed, the result is the Gram completion of the symbolic van Vleck\n  expansion with its positivity conditions; otherwise the error names the order that\n  needs numeric values.\n- The finite generator and its Kossakowski matrices are assembled in the exact native\n  algebra and lifted term by term, so large exact rationals no longer overflow SQA's\n  Rational{Int} in dissipator products.\n- The BigInt fallback also catches InexactError from Int128 conversion.\n- Regression tests: Gram completion where S = 0, drive-phase independence including the\n  newborn channel, a dense Floquet-propagator oracle for a three-level ladder, and the\n  symbolic route.\n\n* feat: terminal virtual gauge resolves the order-six Kerr obstruction\n\nAt the last retained order the static gauge never has to be stored. For one bosonic site with a\nnumber-diagonal averaged generator plus loss and no resonances, the neutral residual is\nsolvable exactly when the lattice invariants l(Y) = <L|Y|R> of ad L0 vanish. The virtual step\nimposes l = 0 on the neutral polynomial parameters, cancels charged entries, passes the\nmixed-charge cross units as gauge images to the exact static solve, and records the order as\nvirtual. The charged ladder of [L0, S] = Y never terminates there, so no polynomial or finite\nrational gauge exists.\n\nDriven Kerr with loss at order six now births a' with rate 18/78125 and a'^2 a^2 with rate\n5832/78125 per omega^6, equal for Bloch-Feshbach and Hori-Deprit. The static part of the\norder-N micromotion is a formal operator and is omitted; the realization lists the virtual\norders. A non-polynomial gauge below the last order still throws.\n\nAlso document the supported algebras and the symbolic route in the GKSLNormalForm docstring,\nand refresh the ratchet baselines for the new files.\n\n* perf: halve native first-call latency and triple steady-state speed\n\nThe first GKSLNormalForm call on driven Kerr at order five drops from about 30 s to about 16 s\nand the steady-state call from 0.9 s to 0.32 s.\n\n- The BigInt retry runs through invokelatest, so inference no longer compiles the whole\n  Rational{BigInt} pipeline when the Int128 attempt succeeds. The expansion is consumed\n  behind a dispatch barrier.\n- Exact matrix kernels are explicit loops instead of generic broadcasting and LinearAlgebra\n  matmul on rational matrices, which also removes most of the runtime cost.\n- The charge nullspace runs over the real scalar type instead of an extra Rational{BigInt}.\n- Tests cover the BigInt fallback and agreement of the loop kernels with the generic ones.\n\n* feat: exact symbolic static gauges by verified reconstruction\n\nWith symbolic parameters and a nonzero or virtual static gauge, GKSLNormalForm now returns the\nexact symbolic expansion instead of throwing. The exact solver is sampled at small rational\npoints in the sign region of a generic reference point. Every sample must reproduce the\nreference structure: frame, channel onsets, static and virtual orders, and nonzero pattern.\nEach output is rebuilt as a rational function over Q(i) from per-variable degree detection,\na homogeneity ray, and an exact solve on the detected support. The result is verified exactly\nat the reference point and two fresh points before it is returned. The lifting code takes a\ncoefficient map, so the symbolic result reuses the numeric assembly.\n\nSymbolic driven Kerr with loss at order five births a^2 - 2a'a with rate 2 F^2 kappa chi^2 per\nomega^4, which equals the exact numeric expansion when substituted, with BF = HD. The positivity\nconditions state the sampled sign region and the channel rates. The faster exact-field route\nvia Nemo is tracked in #389.\n\nAlso replace a closure that captured a @nospecialize argument in algebraic_gauge_directions\nwith a parametric callable struct, which fixes the CheckConcreteStructs failure on #390.\n\n* fix: detect structurally negative symbolic leading rates as obstructions\n\nstructural_sign now classifies a negative exact rational times even powers of real symbols and asserted-nonnegative factors as nonpositive. A negative retained onset therefore raises CompletionObstruction in both Gram and Spectral, instead of being recorded as a positivity condition (Spectral) or misreported as a Puiseux onset (Gram). Closes #379.\n\n* refactor: extract asserted_condition_sign to keep structural_sign under the complexity ratchet\n\n* test: cover each branch of the structural negative-product sign rule\n\n* fix: Spectral completion without square roots or gcd simplification\n\nBranch vectors are kept unnormalized with their squared norm stored separately, so completed coefficients contain no square roots of rational functions. Matrix entries are sums of per-branch quotients with expanded numerators, Hermiticity is first checked structurally, and floating roundoff on symbolic monomials is recognised in GKSL coordinate extraction. Closes #380.\n\n* refactor: annotate and loop roundoff helpers to satisfy the complexity ratchet\n\n* fix: canonicalise completion scalars to trigonometric form\n\ncompletion_scalar now applies trigonometric_form before to_num, so conj(exp(iφ)) = exp(-iφ) for real φ is recognised and Hermitian series with symbolic phases complete. Closes #381.\n\n* docs: state how the dissipative frame selects the positive continuation\n\nGram is the perturbative Cholesky gauge: invariant under lower-triangular frame changes, including rescaling, but dependent on frame ordering. Spectral depends on the coordinate inner product. Addresses item 1 of #382.\n\n* build: track SecondQuantizedAlgebra 0.13 and drop the exact-coefficient workaround\n\nSQA 0.13 keeps coefficients exact by type, so Liouvillian no longer promotes\nGaussian-integer factors before they meet a raw symbolic coefficient. The\nnative lowering reads the new native slot and Gaussian-rational scalars.\n\n* refactor: remove the unreachable AlgebraicChargeInverse route; test phase-site models\n\nThe driver solves charged static equations with native_exact_static_solve, so\nAlgebraicChargeInverse, sparse_exact_solve and their helpers were reached only\nby one internal test. Position and momentum models were not exercised at all\nand now have a public-API test.\n\n* refactor: slim GKSLNormalForm to the exact native core\n\nDrop the research scaffolding that production never reached (dense, exact and\ncharge-graded representations, the float static solve, the shared Gram kernel,\nthe HomologicalInverse seam, test/gksl prototypes) and defer symbolic\nparameters and the terminal virtual gauge to follow-up work.\n\nReject what cannot be expanded exactly instead of continuing silently: floats,\nsymbolic parameters, operators with different names on one Hilbert space, a\nnon-GKSL averaged generator, order < 1, and an order needing a\nnon-polynomial static gauge.\n\nStrip added code comments, trim the manual page and docstring, and update\nADR 0011, CONTEXT.md and the architecture table.\n\n* refactor: update terms function to use Symbolics namespace and improve documentation\n\n* no nothing\n\n* refactor: improve exact numeric radical functions and add condition_sign utility\n\n* more proper docs\n\n* code review\n\n* inference\n\n* inference",
          "timestamp": "2026-10-11T18:16:33+02:00",
          "tree_id": "f2c032a0abf1d89df88f8f8ff124fa58744ec1df",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/c4ab31753a288ba746f11732412fe1db4604f5b8"
        },
        "date": 1791735726578,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1421953,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2541040\nallocs=21406\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 5539074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13573752\nallocs=95773\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 3738219.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9945504\nallocs=70080\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 56278992,
            "unit": "ns",
            "extra": "gctime=11471592\nmemory=94586656\nallocs=635786\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 471250.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329272\nallocs=5435\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 962816,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1033792\nallocs=13268\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 714882,
            "unit": "ns",
            "extra": "gctime=0\nmemory=584936\nallocs=9220\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1812732.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2507632\nallocs=36838\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1038757,
            "unit": "ns",
            "extra": "gctime=0\nmemory=815792\nallocs=13224\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 6972515.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8356192\nallocs=97556\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2361027,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2437624\nallocs=38124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 40621870,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43814616\nallocs=581925\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 386071.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=199440\nallocs=4000\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 488456,
            "unit": "ns",
            "extra": "gctime=0\nmemory=344416\nallocs=5991\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 663320,
            "unit": "ns",
            "extra": "gctime=0\nmemory=629400\nallocs=9904\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 768422.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=390208\nallocs=8264\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1076412,
            "unit": "ns",
            "extra": "gctime=0\nmemory=812400\nallocs=14255\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2434444,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2507288\nallocs=39166\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 369887.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158528\nallocs=3662\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 744126.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330448\nallocs=7774\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 7561226.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3049944\nallocs=71518\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 6228417,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2798888\nallocs=61329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1681795,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1791328\nallocs=18557\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1592830,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1483904\nallocs=16952\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 11856724,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4738264\nallocs=118634\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "15a8f71f9f48186a76f2e15f7c857a8f75161e44",
          "message": "feat: mark completion exceptions as public (#336)\n\nCompletionObstruction and FractionalJumpOnset are thrown by\npositive_completion, but callers could only catch them through unmarked\ninternal names. Mark both @public so they have a stable qualified name, and\nimport them into the test workers together with the rest of the expert API,\nsince the @public list and that import list change together.\n\nBoth docstrings now name the throwing function and describe their fields,\nand FractionalJumpOnset notes that Spectral() can represent an odd onset.\nThe completion manual lists both under a new \"When completion fails\"\nsection.",
          "timestamp": "2026-10-11T21:52:43+02:00",
          "tree_id": "15b4023e572350e80380057b4825b683a6810e31",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/15a8f71f9f48186a76f2e15f7c857a8f75161e44"
        },
        "date": 1791748706402,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1859689.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2541040\nallocs=21406\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 6407794.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13573752\nallocs=95773\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 4539539,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9945504\nallocs=70080\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 61186806,
            "unit": "ns",
            "extra": "gctime=12474478.5\nmemory=94586656\nallocs=635786\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 672822,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329272\nallocs=5435\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1240298,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1033792\nallocs=13268\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 834085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=584936\nallocs=9220\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 2119117,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2507632\nallocs=36838\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1449460,
            "unit": "ns",
            "extra": "gctime=0\nmemory=815792\nallocs=13224\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 8218179.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8356032\nallocs=97556\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2917856,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2437624\nallocs=38124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 45464385,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43814872\nallocs=581925\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 589661,
            "unit": "ns",
            "extra": "gctime=0\nmemory=199440\nallocs=4000\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 705739,
            "unit": "ns",
            "extra": "gctime=0\nmemory=344416\nallocs=5991\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 896322,
            "unit": "ns",
            "extra": "gctime=0\nmemory=629400\nallocs=9904\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 1154376,
            "unit": "ns",
            "extra": "gctime=0\nmemory=390208\nallocs=8264\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1510074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=812400\nallocs=14255\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 3035595,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2507288\nallocs=39166\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 569659,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158528\nallocs=3662\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 1120603,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330448\nallocs=7774\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 10834963.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3049944\nallocs=71518\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 8868238.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2798888\nallocs=61329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 2246406,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1791328\nallocs=18557\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2138915,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1483904\nallocs=16952\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 16893326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4738264\nallocs=118634\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "9111730920f0a33ff184e36e990fdeea0b93c256",
          "message": "feat: exact symbolic parameters in GKSLNormalForm (Nemo) (#397)\n\n* feat: exact symbolic parameters in GKSLNormalForm over a Nemo rational-function field\n\nGKSLNormalForm now expands generators whose coefficients are polynomials in up\nto eight symbolic parameters. The native recurrence runs over Q(p1, ..., pk)\nwith Nemo's canonical fractions, so zero tests are exact. Parameters are\npositive unless listed in the new negative keyword; a sign the solver decides\nmust be fixed throughout that region, otherwise an ArgumentError names it.\npositivity_conditions lists the signed parameters, the drive frequency and the\nchannel rates. Nemo becomes a direct dependency; ADR 0012 records why.\n\n* fix: review findings on symbolic GKSLNormalForm; drop eager zero construction\n\n- reject the drive frequency inside a coefficient instead of treating it as\n  an order-one parameter\n- fall through to the next static-gauge family on an undecidable sign, as on\n  a missing PSD lift\n- store negative as a tuple so equal gauges compare equal\n- regularity_conditions returns the irreducible denominators of the returned\n  coefficients; positivity_conditions holds only the declared signs\n- hash a constant ParameterScalar like the equal Rational\n- build the ring with the public polynomial_ring\n- simplify the scalar wrapper: real fast paths in complex inv and division,\n  drop unused Base methods, share the polynomial-coefficient check\n- accumulate!, native_gauge_metric and hermitian_part no longer build a zero\n  or convert a constant per entry\n\n* style: drop unused argument names in coefficient-domain dispatch\n\n* perf: lift each channel operator with one simplify instead of one per order",
          "timestamp": "2026-10-11T21:55:54+02:00",
          "tree_id": "3735b2c91332592f45597171ce9ab449a42c6a90",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/9111730920f0a33ff184e36e990fdeea0b93c256"
        },
        "date": 1791749017230,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1819261,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2542096\nallocs=21394\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 6213504,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13574808\nallocs=95761\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 4416776,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9946560\nallocs=70068\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 59810495,
            "unit": "ns",
            "extra": "gctime=11108615\nmemory=94587712\nallocs=635774\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 671499,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329400\nallocs=5427\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1229097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1033920\nallocs=13260\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 829601,
            "unit": "ns",
            "extra": "gctime=0\nmemory=585064\nallocs=9212\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 2103481,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2507760\nallocs=36830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1440132,
            "unit": "ns",
            "extra": "gctime=0\nmemory=815984\nallocs=13208\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 8313916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8356000\nallocs=97540\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2833141.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2437816\nallocs=38108\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 44367057.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43814920\nallocs=581909\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 584282,
            "unit": "ns",
            "extra": "gctime=0\nmemory=199568\nallocs=3992\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 698846,
            "unit": "ns",
            "extra": "gctime=0\nmemory=344544\nallocs=5981\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 886146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=629560\nallocs=9892\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 1141584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=390400\nallocs=8248\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1484640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=812592\nallocs=14237\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2953631,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2507512\nallocs=39146\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 563694,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158656\nallocs=3654\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 1099114.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330640\nallocs=7758\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 10341056,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3039096\nallocs=71248\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 8436381,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2748360\nallocs=59956\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 2175932,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1797648\nallocs=18739\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 2084763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1490224\nallocs=17134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 16238525,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4760088\nallocs=119362\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "d24a6f22064db5909f8cd533a2c550209b5579c8",
          "message": "feat: Floquet–Kraus period map (#398)\n\n* feat: Floquet-Kraus period map\n\nAdd period_map, which expands the one-period map of a driven Lindblad\nsystem as a completely positive superoperator built from exact jump\namplitudes in the van Vleck frame of the no-jump generator. The native\nBloch-Feshbach and Hori-Deprit recurrence now takes a pluggable static\nslot so that the period map reuses it with a zero van Vleck slot.\n\n* feat: continuous Kraus form and covariance of the period map\n\nExpose the retained jump amplitudes as a continuous Kraus family with\njump_amplitudes, and state and test the covariance of the map under\nchannel mixing and affine Lindblad shifts.\n\n* fix: concrete types at the SQA coefficient boundary\n\nSQA.to_num and isone(::QAdd) now infer abstract return types, which\nbroke the completion optimizer-stability JET checks on main. Assert\nComplex{Num} and Bool where the package reads them.\n\n* Refactor and optimize code across multiple files\n\n- Moved the `inverse_drive_power` function from `conversion.jl` to `periodic_generator.jl` for better organization.\n- Introduced `kick_inverse` function in `no_jump_frame.jl` to streamline kick operations.\n- Enhanced `lift_graded_superoperator` to accept precomputed scales, improving performance.\n- Updated `period_map` to utilize `JumpRecords` for better management of jump-related data.\n- Replaced dense operator construction with a shared `tomatrix` function in various test files for consistency and reduced code duplication.\n- Improved error handling in `checked_native_parameters` to enforce limits on symbolic parameters.\n- Refactored `floquet_recurrence` functions to accept algorithm as a parameter, enhancing flexibility.\n- Simplified kernel operations in `time_kernel.jl` by consolidating monomial handling functions.\n- Enhanced test cases to utilize new shared functions for operator matrix construction.\n\n* feat: enhance record_map! and record! functions for improved handling of diagonal and upper matrices\n\n* refactor: period map record targets, shared sandwich helpers and clearer names\n\nGroup the record accumulators in RecordTargets, split the off-diagonal Gram row out of record_map!, and route scale_operator and sandwich through the in-place helpers. Add a scaled add_sandwich! for the diagonal Gram term. Rename graded_value to endpoint_value and open_upper to upper_is_jump. Output is unchanged.\n\n* fix: SQA 0.14 compat and a non-deprecated Nemo ring constructor\n\nSQA main is now 0.14, so widen the compat to 0.13, 0.14. Nemo 0.57 deprecates polynomial_ring_only, and the warning broke the doctests on Julia 1.10; use AbstractAlgebra's poly_ring, which exists without deprecation under both Nemo 0.56 and 0.57.\n\n* fix: keep SQA 0.13 in compat for Julia 1.10\n\nJulia 1.10 ignores [sources] and resolves the registered SQA, whose latest release is 0.13.\n\n* better API\n\n* docs\n\n* fix: qadd returns a QAdd under SQA operator function expressions\n\nSQA main now has formal operator function expressions (QExpr <: QField), so 1 * x on a QField infers Union{QAdd,QExpr} and JET flags every caller. qadd now checks that the product is a QAdd and rejects other expressions with an ArgumentError.\n\n* fix: type the SQA boundaries by the supported operator types\n\nSQA main adds formal operator function expressions, QExpr <: QField, which the package does not support. Signatures now take OperatorExpression = Union{QSym,QAdd} instead of QField, and qadd has one method per supported type, so unsupported expressions fail by dispatch instead of a runtime check. SQA #294 makes to_num and isone(::QAdd) infer concrete types, so the coefficient_number wrapper and the Bool assertion of 6b755de are removed.",
          "timestamp": "2026-10-12T11:53:45+02:00",
          "tree_id": "fa891fbbd400854a022161915f62c9cd7c172fea",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/d24a6f22064db5909f8cd533a2c550209b5579c8"
        },
        "date": 1791799202318,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1268127,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2540080\nallocs=21327\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 5948904,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13571096\nallocs=95694\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 3197963,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9943968\nallocs=70017\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 48843341,
            "unit": "ns",
            "extra": "gctime=11994481\nmemory=94524160\nallocs=635099\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 388654.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=327448\nallocs=5369\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 762977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1030144\nallocs=13164\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 516366,
            "unit": "ns",
            "extra": "gctime=0\nmemory=579752\nallocs=9084\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1658749,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2476048\nallocs=36152\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 764464.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=810224\nallocs=13048\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 5589282,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8250544\nallocs=95225\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 1686650.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2422264\nallocs=37532\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 29284276.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42569296\nallocs=551775\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 287637,
            "unit": "ns",
            "extra": "gctime=0\nmemory=198000\nallocs=3942\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 349510.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=341824\nallocs=5907\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 510738,
            "unit": "ns",
            "extra": "gctime=0\nmemory=623576\nallocs=9750\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 643382.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=388480\nallocs=8172\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 796558,
            "unit": "ns",
            "extra": "gctime=0\nmemory=801552\nallocs=13971\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 1668037.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2490952\nallocs=38549\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 275940,
            "unit": "ns",
            "extra": "gctime=0\nmemory=157376\nallocs=3610\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 557553.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329200\nallocs=7692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 5225396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3064088\nallocs=71851\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 4136508.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2804280\nallocs=61361\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1135738,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1793680\nallocs=18570\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1071330.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1487280\nallocs=17013\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 8015500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4779912\nallocs=119770\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "313a45cfb650d55c1d3cabdb61456e6792e2534e",
          "message": "build: drop the numeric-radical workaround now that SQA keeps radicals exact (#407)",
          "timestamp": "2026-10-12T15:07:17+02:00",
          "tree_id": "fde8cd93c8ff516c16227f47d1ade2c3c1116deb",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/313a45cfb650d55c1d3cabdb61456e6792e2534e"
        },
        "date": 1791810829909,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1430679,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2540080\nallocs=21327\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 5416755.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13573432\nallocs=95686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 3676115.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9946272\nallocs=70033\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 58413902.5,
            "unit": "ns",
            "extra": "gctime=13068224\nmemory=94566752\nallocs=635307\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 455247,
            "unit": "ns",
            "extra": "gctime=0\nmemory=327448\nallocs=5369\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 919878,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1030144\nallocs=13164\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 690661,
            "unit": "ns",
            "extra": "gctime=0\nmemory=579752\nallocs=9084\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1638046,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2476048\nallocs=36152\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1027279,
            "unit": "ns",
            "extra": "gctime=0\nmemory=809744\nallocs=13040\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 6759001,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8243056\nallocs=95189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2107224,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2410552\nallocs=37340\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 35522722.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42719016\nallocs=553411\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 392331.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=198000\nallocs=3942\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 475418,
            "unit": "ns",
            "extra": "gctime=0\nmemory=341824\nallocs=5907\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 649886,
            "unit": "ns",
            "extra": "gctime=0\nmemory=623576\nallocs=9750\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 782981,
            "unit": "ns",
            "extra": "gctime=0\nmemory=388480\nallocs=8172\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1042573,
            "unit": "ns",
            "extra": "gctime=0\nmemory=803472\nallocs=14003\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2228792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2479240\nallocs=38357\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 380909.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=157376\nallocs=3610\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 760943,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329200\nallocs=7692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 6795241,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2920984\nallocs=68045\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 5348688,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2533096\nallocs=54562\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1465385.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1750360\nallocs=17439\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1395219.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1443960\nallocs=15882\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 10924871,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4615544\nallocs=115404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "2f3b4212eeb1c4223a9636e54d3601b658c5696d",
          "message": "fix: read SQA monomial scalars through monomial_scalar (#405)\n\n* fix: read SQA monomial scalars through monomial_scalar\n\n* fix: read the monomial scalar with SQA.monomial_scalar and require SQA 0.14",
          "timestamp": "2026-10-12T16:54:07+02:00",
          "tree_id": "1e0d282cd1bd44cb826a8fe1c19fe86b6f5c3eac",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/2f3b4212eeb1c4223a9636e54d3601b658c5696d"
        },
        "date": 1791817527217,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1524103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2567504\nallocs=21327\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 5879333.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13657800\nallocs=95734\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 3885768,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10043008\nallocs=70049\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 59451103,
            "unit": "ns",
            "extra": "gctime=16930783.5\nmemory=94916336\nallocs=635355\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 458691,
            "unit": "ns",
            "extra": "gctime=0\nmemory=327512\nallocs=5369\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1116851,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1031936\nallocs=13164\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 664592.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=581032\nallocs=9084\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1851449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2476560\nallocs=36152\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1085123,
            "unit": "ns",
            "extra": "gctime=0\nmemory=807824\nallocs=13040\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 7783035.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8220256\nallocs=95189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2338249,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2385912\nallocs=37340\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 37804942,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41984672\nallocs=551669\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 389998,
            "unit": "ns",
            "extra": "gctime=0\nmemory=197424\nallocs=3942\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 481965,
            "unit": "ns",
            "extra": "gctime=0\nmemory=342080\nallocs=5907\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 663174,
            "unit": "ns",
            "extra": "gctime=0\nmemory=624952\nallocs=9750\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 786800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=387392\nallocs=8172\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1090441,
            "unit": "ns",
            "extra": "gctime=0\nmemory=797072\nallocs=14003\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2468164,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2454472\nallocs=38357\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 371945,
            "unit": "ns",
            "extra": "gctime=0\nmemory=157184\nallocs=3610\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 770373.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328752\nallocs=7692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 7636803,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2954776\nallocs=68855\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 6085222.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2598888\nallocs=56182\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1701919,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1750024\nallocs=17427\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1596569,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1443624\nallocs=15870\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 12240216.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4649336\nallocs=116214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "ef9e588f165f5dcfe995bf08a4687e481950ca60",
          "message": "Update Julia version in CI workflow (#416)",
          "timestamp": "2026-10-12T21:07:52+02:00",
          "tree_id": "ee1c08f9c291613c35e21d8ee17bfc48110a3583",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/ef9e588f165f5dcfe995bf08a4687e481950ca60"
        },
        "date": 1791832776771,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1822116,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2567504\nallocs=21327\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 6034338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13654104\nallocs=95678\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 4333121.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10036640\nallocs=69937\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 57325486,
            "unit": "ns",
            "extra": "gctime=11042754.5\nmemory=94880368\nallocs=635019\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 673854.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=327512\nallocs=5369\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1216308.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1031936\nallocs=13164\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 801263,
            "unit": "ns",
            "extra": "gctime=0\nmemory=581032\nallocs=9084\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1948904,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2476560\nallocs=36152\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1437723,
            "unit": "ns",
            "extra": "gctime=0\nmemory=807824\nallocs=13040\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 7886797,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8221616\nallocs=95189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2726296,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2379256\nallocs=37212\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 38417642.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41796456\nallocs=548318\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 584327.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=197424\nallocs=3942\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 691262,
            "unit": "ns",
            "extra": "gctime=0\nmemory=342080\nallocs=5907\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 870903.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=624952\nallocs=9750\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 1184845,
            "unit": "ns",
            "extra": "gctime=0\nmemory=387392\nallocs=8172\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1475448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=797072\nallocs=14003\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2887838.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2447816\nallocs=38229\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 563553,
            "unit": "ns",
            "extra": "gctime=0\nmemory=157184\nallocs=3610\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 1109694.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328752\nallocs=7692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 9880790,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2923576\nallocs=68045\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 7680981.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2536488\nallocs=54562\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1979952,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1750024\nallocs=17427\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1901004,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1443624\nallocs=15870\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 15733760.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4618136\nallocs=115404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "2187b77f08f21247839ecc471557bf4639021857",
          "message": "ci: require Julia 1.13 and test on the latest 1.x (#419)\n\n* ci: require Julia 1.11 and test it as the lower bound\n\n* ci: require Julia 1.13 and test only 1.13\n\n* ci: test on the latest Julia 1.x",
          "timestamp": "2026-10-13T09:10:59+02:00",
          "tree_id": "4ec049624f3da0694d46f8879736d6b2b69c3451",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/2187b77f08f21247839ecc471557bf4639021857"
        },
        "date": 1791876120705,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1659709,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2567504\nallocs=21327\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 6066745,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13649592\nallocs=95614\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 4088696,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10036160\nallocs=69937\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 57089518,
            "unit": "ns",
            "extra": "gctime=13189090\nmemory=94871696\nallocs=635019\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 560396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=327512\nallocs=5369\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1181753,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1031936\nallocs=13164\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 744879,
            "unit": "ns",
            "extra": "gctime=0\nmemory=581032\nallocs=9084\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1932880,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2476560\nallocs=36152\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1248277,
            "unit": "ns",
            "extra": "gctime=0\nmemory=807824\nallocs=13040\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 7876205,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8224096\nallocs=95181\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2548392,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2380920\nallocs=37244\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 39260845,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41812120\nallocs=549283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 461698,
            "unit": "ns",
            "extra": "gctime=0\nmemory=197424\nallocs=3942\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 591848,
            "unit": "ns",
            "extra": "gctime=0\nmemory=342080\nallocs=5907\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 789638,
            "unit": "ns",
            "extra": "gctime=0\nmemory=624952\nallocs=9750\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 951593.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=387392\nallocs=8172\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1274837,
            "unit": "ns",
            "extra": "gctime=0\nmemory=795408\nallocs=13971\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2653407,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2449480\nallocs=38261\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 440408,
            "unit": "ns",
            "extra": "gctime=0\nmemory=157184\nallocs=3610\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 924835,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328752\nallocs=7692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 8691853,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2954776\nallocs=68855\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 7010451,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2598888\nallocs=56182\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1954507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1750024\nallocs=17427\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1875347,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1443624\nallocs=15870\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 14005033.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4649336\nallocs=116214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "0b7ada441909855ce539fa69a5c2cb11c1091773",
          "message": "refactor: remove the zero_component field of PeriodicGenerator (#413)\n\n* refactor: remove the zero_component field of PeriodicGenerator\n\n* test: drop the local path source from test/Project.toml",
          "timestamp": "2026-10-13T09:20:17+02:00",
          "tree_id": "ded2ee8af47529e831a0136dc7a4c5849c016319",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/0b7ada441909855ce539fa69a5c2cb11c1091773"
        },
        "date": 1791876382662,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1390826.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2579488\nallocs=21413\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 5447404.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13667032\nallocs=95776\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 3702057,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10047312\nallocs=70008\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 53081535,
            "unit": "ns",
            "extra": "gctime=10911403\nmemory=94889072\nallocs=635069\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 457017.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328056\nallocs=5373\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 924531,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1032448\nallocs=13168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 577450,
            "unit": "ns",
            "extra": "gctime=0\nmemory=575976\nallocs=8910\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1520407,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2458064\nallocs=35552\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1018368,
            "unit": "ns",
            "extra": "gctime=0\nmemory=808368\nallocs=13044\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 6713001,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8221840\nallocs=95193\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2096194,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2374200\nallocs=37038\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 32886717,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41778376\nallocs=547718\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 395342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=197968\nallocs=3946\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 472290,
            "unit": "ns",
            "extra": "gctime=0\nmemory=339888\nallocs=5833\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 613207.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=617736\nallocs=9526\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 780368.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=387936\nallocs=8176\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1027011,
            "unit": "ns",
            "extra": "gctime=0\nmemory=794880\nallocs=13929\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2178475,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2440600\nallocs=38005\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 371316.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=157760\nallocs=3614\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 755943,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329328\nallocs=7696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 6952398,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2954776\nallocs=68855\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 5597305,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2598888\nallocs=56182\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1473623,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1750008\nallocs=17427\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1393860,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1443624\nallocs=15870\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 11253060,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4649336\nallocs=116214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "0b9a73962eeabda6d174cf9a6f27d3e7e0d6bec7",
          "message": "perf: grade period-map time integrals by powers of the period instead of a π field variable (#415)\n\n* perf: grade period-map time integrals by powers of the period instead of a π field variable\n\n* perf: keep parameter period maps out of the Int128 overflow retry\n\nExact maps try Complex{Rational{Int128}} and fall back to BigInt; maps with\nsymbolic parameters dispatch straight to ParameterScalar, so a failure there is\nno longer computed twice.\n\nTests: eight symbolic parameters expand and a ninth is rejected; a rate beyond\nInt128 goes through the BigInt fallback and matches the nearby exact map.",
          "timestamp": "2026-10-13T09:44:06+02:00",
          "tree_id": "f39fe8b25c89b1daa7c3d25365fcf76045343e45",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/0b9a73962eeabda6d174cf9a6f27d3e7e0d6bec7"
        },
        "date": 1791878105662,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1485732,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2579920\nallocs=21434\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 5475561,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13670184\nallocs=95845\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 3704080,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10054688\nallocs=70157\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 55732070,
            "unit": "ns",
            "extra": "gctime=13225585\nmemory=94926848\nallocs=635554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 450631.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328344\nallocs=5387\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1100920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1032736\nallocs=13182\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 589884,
            "unit": "ns",
            "extra": "gctime=0\nmemory=576264\nallocs=8924\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1704301,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2458352\nallocs=35566\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1071540.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=808528\nallocs=13064\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 7419419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8218640\nallocs=95153\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2247186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2347960\nallocs=36554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 36672672.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41565136\nallocs=545055\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 381487,
            "unit": "ns",
            "extra": "gctime=0\nmemory=198256\nallocs=3960\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 470291,
            "unit": "ns",
            "extra": "gctime=0\nmemory=340240\nallocs=5849\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 642236,
            "unit": "ns",
            "extra": "gctime=0\nmemory=618152\nallocs=9544\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 771566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=388512\nallocs=8204\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1069067,
            "unit": "ns",
            "extra": "gctime=0\nmemory=793856\nallocs=13927\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2358087,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2414488\nallocs=37525\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 366462,
            "unit": "ns",
            "extra": "gctime=0\nmemory=158048\nallocs=3628\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 754497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=329904\nallocs=7724\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 7364146.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2963048\nallocs=69121\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 5987524.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2606568\nallocs=56444\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1696575.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1751672\nallocs=17479\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1597910,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1445288\nallocs=15922\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 11803807.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4662104\nallocs=116621\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
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
          "id": "e1ba0f0501f050501ccae5115d850f43e568b6f3",
          "message": "fix: drop use of removed SecondQuantizedAlgebra.is_constant_term (#421)\n\n* fix: test constant Poly terms through Monomial fields after SecondQuantizedAlgebra dropped is_constant_term\n\n* chore: restore SecondQuantizedAlgebra source",
          "timestamp": "2026-10-16T14:20:15+02:00",
          "tree_id": "06ae3143920d571ae3ae0f346bd978759664616a",
          "url": "https://github.com/oameye/FloquetExpansions.jl/commit/e1ba0f0501f050501ccae5115d850f43e568b6f3"
        },
        "date": 1792153908222,
        "tool": "julia",
        "benches": [
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 3",
            "value": 1623518,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2535888\nallocs=20333\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Bloch-Feshbach order 5",
            "value": 5795762,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13604088\nallocs=93872\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 3",
            "value": 4088878.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10012272\nallocs=69140\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven dissipative qubit/Hori-Deprit order 5",
            "value": 53873687.5,
            "unit": "ns",
            "extra": "gctime=9754314\nmemory=94902112\nallocs=634721\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 3",
            "value": 552310,
            "unit": "ns",
            "extra": "gctime=0\nmemory=299936\nallocs=4694\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Bloch-Feshbach order 5",
            "value": 1071107,
            "unit": "ns",
            "extra": "gctime=0\nmemory=981016\nallocs=11635\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 3",
            "value": 664280,
            "unit": "ns",
            "extra": "gctime=0\nmemory=540560\nallocs=7964\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Driven qubit/Hori-Deprit order 5",
            "value": 1748051,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2364096\nallocs=32466\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 3",
            "value": 1226201.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=741888\nallocs=11358\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Bloch-Feshbach order 5",
            "value": 7621042.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7917328\nallocs=84931\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 3",
            "value": 2518578.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2268768\nallocs=34134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Expansion Algorithm/Kerr parametric oscillator/Hori-Deprit order 5",
            "value": 39996079,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41006312\nallocs=523101\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 1",
            "value": 483271,
            "unit": "ns",
            "extra": "gctime=0\nmemory=173360\nallocs=3396\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 2",
            "value": 577056.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=312400\nallocs=5178\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Driven qubit/order 3",
            "value": 735597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=582128\nallocs=8574\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 1",
            "value": 958807,
            "unit": "ns",
            "extra": "gctime=0\nmemory=336512\nallocs=7040\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 2",
            "value": 1253417.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=733024\nallocs=12456\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Floquet Expansion/Kerr parametric oscillator/order 3",
            "value": 2580731,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2334848\nallocs=35091\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Driven qubit/symbolic input",
            "value": 464526.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=133152\nallocs=3064\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Fourier Expansion/Kerr parametric oscillator/symbolic input",
            "value": 934172,
            "unit": "ns",
            "extra": "gctime=0\nmemory=277904\nallocs=6560\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Gram fixed frame",
            "value": 9789672.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2923832\nallocs=67876\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Driven qubit/Spectral fixed frame",
            "value": 7495202,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2519920\nallocs=54066\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram automatic frame",
            "value": 1982157.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1749912\nallocs=17424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Full-rank bosonic/Gram fixed frame",
            "value": 1894404,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1443528\nallocs=15867\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Positive Completion/Recursive dark onset/Gram fixed frame",
            "value": 15696980,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4618184\nallocs=115231\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}