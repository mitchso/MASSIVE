# Data requirements
Data for each analysis should be provided as a folder of text files, where each text file corresponds to an individual mass spectrum.

Each text file should be named according to the following convention to ensure proper parsing:
     [year]\_[month]\_[day]\_[run#]\_[chip#]\_[well]\_[replicate].txt \
i.e. 2025_12_06_0001_0_E1_1.txt

Where **[chip#]** refers to the 'chip' on a Bruker MALDI target (0 for regular positions, 1 for calibrant spots) and **[well]** refers to the position on the target (i.e. A1, G12, L23, etc). 

The data inside the text file should be two columns without headers (m/z and intensity), separated by a single space:

i.e.
```aiignore
...
3199.09 61672
3199.43 64551
3199.77 67041
3200.11 69142
3200.46 70855
3200.80 72179
3201.14 73115
3201.48 73662
...
```

See /data for example files.