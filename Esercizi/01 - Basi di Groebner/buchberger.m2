sCoppia = (f,g) -> (
	lf := leadTerm(f);
	lg := leadTerm(g);
	xg := lcm(leadMonomial(f), leadMonomial(g));
	xg//lf * f - xg//lg * g
)

divisione = (f, divisori) -> (
	quozienti := new MutableList from (divisori-divisori);
	resto := 0;
	restoAux := f;
	while restoAux != 0 do (
		j := position(divisori, (g) -> (leadTerm(restoAux) % ideal(leadTerm(g)) == 0));
		if instance(j,Nothing) then (
			resto = resto + leadTerm(restoAux);
			restoAux = restoAux - leadTerm(restoAux);
		) else (
			quozienti#j = quozienti#j + leadTerm(restoAux)//leadTerm(divisori#j);
			restoAux = restoAux - divisori#j*leadTerm(restoAux)//leadTerm(divisori#j);
		);
	);
	(toList(quozienti), resto)
)

buchberger = (generatori) -> (
	coppie := subsets(generatori, 2);
	nuoviElementi := apply(coppie, x -> sCoppia(x#0,x#1));
	nuoviElementi = unique(nuoviElementi);
	nuoviElementi = select(nuoviElementi, x -> (divisione(x,generatori))#1 != 0);
	if #nuoviElementi == 0 then (
		generatori
	) else (
		buchberger(generatori | nuoviElementi)
	)
)