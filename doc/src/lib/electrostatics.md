# Electrostatics

## Potential types
```@docs
    PotentialType
    SingleLayer
    DoubleLayer
```

## Locality assumption
```@docs
    LocalityType
    LocalES
    NonlocalES
```

## Potentials

### Molecular potentials
```@docs
    φmol
    ∂ₙφmol
    ProteinES.∇φmol
```

### Interior potentials
```@docs
    ProteinES.BEM.φΩ
    ProteinES.Born.φΩ
```


### Exterior potentials
```@docs
    ProteinES.BEM.φΣ
    ProteinES.Born.φΣ
```

## Energies
```@docs
    ProteinES.BEM.rfenergy
```
