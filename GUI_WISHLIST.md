GUI test 14/10/19

# THINGS TO FIX

### schermata *select a file*
- [ ] BUG: se clicco sul link "README.md" nella finestra di Help mi dà un errore perché non trova il file:
```file:///home/lercole/software/README.md: Error when getting information for file “/home/lercole/software/README.md”: No such file or directory```

- [ ] il cambiamento delle settings (lingua o font) avviene solo al riavvio del programma

- [ ] La schermata di selezione dei file confonde un po' le idee: sopra vengono elencati i file presenti nella cartella corrente, sotto c'è il preview del file selezionato, e in fondo nella casella "selected" è scritto tutto il percorso del file aperto. Il nome di questo file è spesso illeggibile perché troppo lungo.
Per rendere tutto più lineare (dall'alto in basso) metterei:

- La casella per la selezione del file all'inizio.
- Non è chiaro perché ci sia una lista di files della cartella corrente, se il file si seleziona tramite menù nella casella precedente.
- Sarebbe bello che ci fosse la possibilità di navigare nelle cartelle direttamente in questa schermata (senza aprire una finestrella - tipo il Finder del Mac), e così si potrebbbe selezionare il file dalla lista di file nella cartella corrente.
- Se non si può: o eliminiamo la lista dei file, oppure la teniamo ma scriviamo qual è la cartella corrente ed *evidenziamo* il file selezionato.

- [ ] Vanno chiariti i tipi di file supportati... Ci dovremmo mettere un link al manuale (quando ci sarà un manuale), ma intanto si può scrivere una descrizione breve dei formati lì.

- [ ] C'è un modo per far sì che se chiudo e riapro il programma vada ad aprire direttamente la cartella dell'ultimo file usato, senza ripartire sempre dalla Home?

### schermata *define headers*
- [x] BUG: Se uso il file binario data/Silica.npy, setto `Temperature: None`, clicco next, poi back: Temperature sparisce :)

- [x] BUG: il programma permette di settare più colonne `Temperature`, `Volume`, o `DT_FS`

- [ ] Vanno chiarite le unità di misura.

### schermata *set variable*
- [ ] permetterei all'utente di cambiare i valori delle variabili *anche* se queste vengono lette da file. Il valore che compare dovrebbe essere quello calcolato

### schermata $f^*$
- [ ] Manca un titolo/spiegazione.

- [ ] Il margine sx del grafico è troppo piccolo, i numeri sull'asse y mi vengono tagliati fuori.

- [ ] a cosa serve il pulsante bianco quadrato a sinistra di Reset View?

- [ ] farei lo sfondo azzurro della regione selezionata dallo slider un po' più chiaro, oppure giallino chiaro

- [ ] il sistema zoom-in reset-view può essere un po' macchinoso. Manca la possibilità di fare zoom-in più volte e zoom-out.

- [x] BUG: se seleziono una zona e clicco zoom-in, poi clicco Back le variabili tornano a zero (se il file è di tipo Table). Oltre a questo, se ora imposto le variabili (o tengo quelle già presenti) e vado avanti, $f^*$ rimane bloccata su zero, e lo slider non funziona. L'unico modo per sistemarlo è impostare manualmente un nuovo $f^*$ nella casella.

- [ ] Se uno decide di utilizzare i bottoni del grafico per zommare/pannare, appena tocca lo slider lo zoom torna a quello di default. Questo comportamento potrebbe essere fastidioso... non so se si può risolvere (ad esempio forse si possono ottenere i nuovi limiti del grafico con `plt.xlim()`).

### schermata $P^*$
- [ ] manca titolo/spiegazione

- [ ] cos'è la linea rossa a $f=0$?

- [ ] IMPORTANT: $P^*$ a differenza di $f^*$ *non* è un parametro che l'utente sceglie, di default. Esso viene calcolato automaticamente (è `aic_Kmin + 1`, mi sembra) tramite l'Akaike Information Criterion.

- Scriverei quindi chiaramente qual è il $P^*$ suggerito dall'AIC.
- Inserirei un bottone per resettarlo a questo valore.
- Se l'utente vuole modificarlo glielo lasciamo fare, ma chiamerei questa casella "$P^*$ correction"

- [ ] prevederei la possibilità di modificare il filtro del plot (Filter width) anche qui, così che uno possa modificarlo per analizzare meglio i risultati.


# GENERAL WISH LIST
- [ ] Vanno inserite delle spiegazioni che aiutino a capire cosa rappresenta ogni schermata, senza aver letto approfonditamente il paper. Ad esempio spiegare come scegliere $f^*$.

- [ ] L'help per come è adesso è difficile da leggere. Avrebbe più senso preparare un help html più dettagliato, con degli screenshot di esempio. Tra le cose che vorrei implementare per la versione 1.0 di st c'è la documentazione su readthedocs, fatta con sphinx.

- [ ] si potrebbe fare che se il file letto ha estensione `.npy` o `.npz` viene riconosciuto come binario (dict) direttamente.

- [ ] Per tutti i plot: inserire un pulsante per visualizzare il log-periodogramma al posto del periodogramma.

- [ ] $P^*$ plot: inserire un pulsante per visualizzare i seguenti plot (in finestre esterne):

- coefficienti cepstrali (`HeatCurrent.plot_ck()`)
- analisi di convergenza kappa(P*)  (`HeatCurrent.plot_kappa_Pstar()` e/o `HeatCurrent.plot_L0_Pstar()`
- analisi di $f^*$ -- più difficile, perché richiede di specificare una lista di $f^*$ o `TSKIP` e può essere lenta. Dovremmo aggiungere delle funzioni che generino delle griglie lineari/logaritmiche, così da rendere la selezione più semplice. Per ora concentriamoci sulle altre cose...

### packaging
- [ ] get rid of `uncertainties` package

- [ ] get rid of `tk-html-widgets` if possible

- [ ] ensure that the code is Python 2 compatible
