1. Izvršiti kontrolu kvaliteta FASTQ fajlova alatom FastQC. Priložiti izvrštaj i diskutovati rezultate.

    **> Paired end 1: Per base sequence content je kriterijum koji je označen kao problematičan i on označava propociju svake baze svih readova. U najboljem slučaju, procenat svih baza bi trebao da bude 25%, ali ovde se taj odnos menja od pozicije 90.**
    
    **> Paired end 2: Per base sequence quality je kriterijum koji je označen kao problematičan i on označava kvalitet baza na svim pozicijama svih readova. To je zato što bazni parovi od pozicije 160 imaju veću devijaciju po pitanju kvaliteta.**

2. Mapirati sekvencirane readove na referentni genom hg38 upotrebom alata BWA.
    
    2a. Koliko je readova uspešno mapirano?

    **> Mapped reads: 3348213**

    2b. Koliko je parova readova mapirano tako da su oba para mapirana?

    **> Mapped pairs: 298607631**

    2c. Nacrtati histogram dužina sekvenciranih fragmenata.
    
    ![](data:image/png;base64,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)

3. Izvršiti obradu dobijenog BAM fajla prema GATK protokolu (markiranje duplikata, rekalibracija kvaliteta baza). Koliki su procenati PCR i optičkih duplikata?

    **> Ukupno duplikata: 480157**
    
    **> PCR duplikati: 100 * 480157 / 6775767 = 7.09%**
    
    **> Optički duplikati: 100 * 0 / 6775767 = 0%**

4. Identifikovati mutacije upotrebom alata Haplotype Caller i filtirtati mutacije predefinisanim filterima prema Broad preporukama.

    4a. Koliko je ukupno mutacija identifikovano, koliko od njih su SNP, a koliko INDEL?

    **> SNP: 14577**
    
    **> INDEL: 1796**

    4b. Koliko mutacija prolazi, a koliko ne prolazi kriterijume filtriranja.

    **> Prolazi: 16368**
    
    **> Ne prolazi: 5**

    4c. Izračunati TiTv odnos pre i posle filtriranja.

    **> Pre: 1.9288728149487644**
    
    **> Posle: 1.9294472361809045**

5. Anotirati mutacije alatom Funcotator. Izbrojati različite vrednosti ClinVar značajnosti.

    **> Benign likely benign: 45**
    
    **> Benign: 247**
    
    **> Not provided: 11**
    
    **> Likely benign: 28**
    
    **> Benign other: 1**
    
    **> Likely pathogenic: 1**
    
    **> Pathogenic: 1**
    
    **> Association: 1**
    
    **> Uncertain significance: 3**
    
    **> Conflicting interpretations of pathogenicity affects: 1**
    
    **> Conflicting interpretations of pathogenicity: 2**
    
    **> Risk factor: 2**
    
    **> Benign likely benign association: 1**
    
6. Svi uzorici sadrže određenu količinu kontaminacije DNK materijalom bakterijskog ili virusnog porekla. Većina ovakvih readova se neće mapirati na ljudski genom. Izvući readove koji nisu mapirani u procesu mapiranja, asemblovati ih alatom abyss, i identifikovati organizam od kojeg potiče najduži skafold upotrebom alata Blast.

    **> Nemapirana sekvenca potiče od bakterije [Bradyrhizobium Japonicum](https://en.wikipedia.org/wiki/Bradyrhizobium_japonicum)**
