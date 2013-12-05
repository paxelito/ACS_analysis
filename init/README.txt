LEGGERE CON ATTENZIONE!!!

++++++++++++++++++++++++++++
PER LANCIARE LE SIMULAZIONI:
++++++++++++++++++++++++++++

1) Caricare il file carness_package.tar.gz sul server nella cartella desiderata

2) Estrarre il contenuto del file

	tar -xvzf carness_package.tar.gz

3) A questo punto trovi due tipi di luncher:

	1- luncher_screen.sh : 

Questo luncher utilizza il comando screen per lanciare i terminali virtuali su cui vengono lanciate le simulazioni. A mio parere se screen è installato sul computer (sudo apt-get install screen in caso per installarlo) è molto più comodo perché permette di recuperare velocemente informazioni sulle simulazioni che stanno andando. 
In calce qualche nota sull’utilizzo degli screen 
 
	2- luncher_background.sh

Questo luncher sfrutta il classico “&” a fine comando per lanciare il comando in background. Benissimo anche questo ma una volta utilizzato per poter controllare lo stato delle simulazioni sei obbligato ad andare nella cartella SIMS e guardare i file log. 


3) Verificare che il file luncher (screen o background) sia eseguibile (755)
	
	ls -l permessi devono essere -rwxr-xr-x
	
   se cosi non fosse eseguire

	chmod 755 luncher_background.sh
		oppure
	chmod 755 luncher_screen.sh

4) Eseguire luncher

	./luncher_screen.sh oppure ./luncher_background.sh

5) A questo punto dovrebbe essere estratta la cartella con le simulazioni, compilato il sorgente e partire le simulazioni su 8 core. 

++++++++++++++++++++++++++++++++++++++++++
PER CONTROLLARE LO STATO DELLE SIMULAZIONI
++++++++++++++++++++++++++++++++++++++++++

Per controllare lo stato delle simulazioni abbiamo 3 modi sempre più approfonditi:
1 - Controllare i processi in corso 
	
	eseguire il comando “top” da terminale e verificare che ci siano 8 processi “carness” in corso

2 - Controllare lo stato dei vari lanciatori (disponibile solo se si è utilizzato screen)

	1. eseguire comando “screen -ls” per visualizzare gli screen attivi 
	2. eseguire comando “screen -r pid” oppure “screen -r nomescreen”. Il pid è il numerino che identifica lo screen
	3. All’interno dello screen è possibile monitorare quale simulazione è in corso
	4. IMPORTANTE!!! per uscire dallo screen digitare “ctrl+a+d”, attenzione a non chiudere in altro modo lo screen altrimenti terminano tutte le simulazioni in corso su quello screen

3 - Controllare i file log delle simulazioni

	entrare nella cartella exps/SIMS/ e digitare il comando 
	“tail -f -n250 CH#_N_condition.log” dove # indica la chimica e N il numero della condizione (1->600, 2->60, 3->6, 4->1). Questo comando mostra in tempo reale l’andamento della simulazione. 
	per uscire dal comando fare semplicemente “ctrl+c”

