## Adding moves
The signature of the `add_move` method is as follows:
``` python
system.add_move(move_name, prob=1.0, box_id=-1, **kwargs)
```
where `box_id=-1` means that the move is assiged to all available boxes. Some examples are:
``` python
system.add_move("trans", prob=0.5, dx=0.1)
system.add_move("avbmc", prob=0.25, particle="Si", r_below=0.9, r_above=3.0, box_id=0)
system.add_move("avbmcmol", prob=0.25, molecule=water_mol, r_below=0.95, r_above=3.0, r_inner=1.3)
```

