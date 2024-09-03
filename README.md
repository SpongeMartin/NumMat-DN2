# Martin Starič

## Gauss-Legendrove kvadrature


Projekt vsebuje funkcije za sestavljeno Gauss-Legendrovo kvadraturno formulo za integrale, ki tečejo od 0 do h. Bolj podrobno je raziskan integral od 0 do 5 sin(x)/x
### Zaganjanje kode

Kodo zaženemo tako da:

1. **Aktiviramo okolje:**
   - Odpri Julia REPL in pojdi v način pkg tako da napišeš `]`.
   - Aktiviraj okolje z ukazom:
     ```julia
     activate DN1
     ```

2. **Uporaba kode:**
   - Sledi primeru v `\test\demo.jl`.

3. **Testi:**
    - Testi so napisani v `\test\runtests.jl`.

### Generiranje .tex datoteke

`.tex` datoteka je generirana s skripto `\docs\makedoc.jl`, ki uporablja paket `Weave.jl`, v tem direktoriju se nahaja tudi poročilo.



