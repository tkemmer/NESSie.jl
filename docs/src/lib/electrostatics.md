# Electrostatics
```@meta
    CurrentModule = NESSie
```

```@index
Pages = ["electrostatics.md"]
```

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
    ∇φmol
```

### Interior potentials
```@docs
    BEM.φΩ
    TestModel.φΩ
```


### Exterior potentials
```@docs
    BEM.φΣ
    TestModel.φΣ
```

## Energies
```@docs
    BEM.rfenergy
```
