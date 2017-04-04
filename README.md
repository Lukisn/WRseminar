WR Seminar
==========
**Modellierung der Populationsbilanz von Aerosolen
mittels der Sektionalmethode**

Overview
--------
```
.
├── README.md
├── __pycache__
├── demo/
│   ├── __init__.py
│   ├── demo_functions.py
│   ├── demo_grid.py
│   ├── demo_method.py
│   └── demo_util.py
├── docs/
├── requirements.txt
├── sectional/
│   ├── __init__.py
│   ├── functions.py
│   ├── grid.py
│   └── methods.py
├── cases/
│   ├── __init__.py
│   ├── case_agg.py
│   ├── case_break.py
│   ├── case_break_agg.py
│   ├── case_gro_agg.py
│   ├── case_growth.py
│   ├── case_growth1.py
│   └── case_growth4.py
├── tests/
│   ├── __init__.py
│   ├── test_grid.py
│   └── test_grid.pyc
└── util/
    ├── __init__.py
    └── cmdline.py
```


TODO
----

* integrate silly demo code into tests
* refactor testcases -> extract commom plotting code


Old Notes
---------

immernoch problematisch:
- "Ränder" in Fixed Pivot
    o Fokus auf Cell Average!
- testcase pure growth?
- no testcase pure nucleation?

Anregungen vom letzten mal:
- Stabilität "courant" zahl?
- nummerische diffusion -> siehe "simple" growth
- Gitterstudie
    o Vergleich äquidistant/geometrisch?
    o verfeinerung geometrisch?