��   �   0     �    
     |     �  �      ,  �   -  7  �  �  7  -      F   .     u     �  7   �  �   �  �   q  �   $  I      �   j  �   �  �   �  �  ~  H   -     v  E   �  �   <    �  >   �  �   (  9   �  �   �  �   �  �      �   �   �   a!  �   �!  l   �"      #     :#     T#     n#     �#     �#     �#     �#     �#  z   $     �$     �$     �$  .   �$  6   %     ;%     N%     b%     g%  !   }%  !   �%  '   �%     �%     	&     )&  *   H&  /   s&  %   �&     �&  /   �&  ,   '     8'  4   N'  $   �'     �'     �'     �'     �'     (      0(      Q(  h   r(  <   �(     )  :   +)  &   f)  $   �)     �)  2   �)      *  $   *  /   B*  I   r*     �*     �*  3   �*  =   +  6   X+  d   �+  [   �+      P,  O   q,  .   �,  /   �,      -  A   ;-  .   }-  0   �-  )   �-     .     .     ).  ;   C.  <   .  8   �.  6   �.     ,/     I/     d/  (   y/  I   �/  !   �/  '   0  '   60  9   ^0     �0     �0  0   �0     �0  <   �0  -   )1  @   W1  8   �1  >   �1  /   2  7   @2  D   x2  5   �2  &   �2  '   3     B3  %   J3     p3     �3  
   �3  
   �3  
   �3  
   �3  
   �3  	   �3  	   �3  	   �3  	   �3  	   �3  	   �3  	   	4  "   4  *   64     a4  A   u4  Q   �4  *   	5  @   45  !   u5     �5    �5  �   �7  �  �8    +=  A   -?  ]   o?     �?     �?  G   �?  �   C@  �   �@  3  �A  R   �B  �   RC  X  �C  �   FE  b  �E  o   OH  �   �H  D   CI  �   �I  �  lJ  N   +L  �   zL  =   mM  �   �M  �   ~N  �   O  *  �O  �   Q  �   �Q  n   �R      FS      gS  !   �S  %   �S  !   �S  "   �S      T  !   6T      XT  z   yT      �T  !   U     7U  2   IU  <   |U  (   �U     �U     �U     V  >    V  6   _V  4   �V  (   �V  -   �V  '   "W  1   JW  F   |W  :   �W     �W  1   X  =   GX     �X  E   �X  8   �X  ,   Y     JY     jY  !   �Y     �Y  *   �Y  *   �Y  �   Z  O   �Z     �Z  `   [  <   n[  -   �[  (   �[  E   \  '   H\  1   p\  8   �\  t   �\     P]     j]  2   �]  M   �]  K   ^  �   O^  |   �^  &   R_  m   y_  4   �_  1   `     N`  B   m`  ?   �`  A   �`  0   2a  
   ca  %   na  '   �a  U   �a  U   b  F   hb  F   �b  2   �b  '   )c     Qc  *   gc  _   �c  #   �c  3   d  3   Jd  S   ~d  -   �d      e  L   e     Pe  I   Ye  �   �e  G   .f  L   vf  O   �f  +   g  9   ?g  L   yg  F   �g  0   h  5   >h  	   th  7   ~h  "   �h     �h  
   �h  
   �h  
   i  
   i  
   i  	   $i  	   .i  	   8i  	   Bi  	   Li  	   Vi  	   `i  &   ji  4   �i     �i  W   �i  a   3j  <   �j  M   �j  &    k     Gk     �   V                       l   /       K       	          4   ^   O   �      y   E   �      k   @   o      p   J           M   h   �       $       G       %   |   ]       R       s   ~              e      N   #          &   �      Q   �          b           
   \      6   c       _           Y       "   i   >   .   S                  w          t           u   D       5   r       X               q   �      )      �           !   1   �   A               (       :   C   F   }   B   H   m   +       '   �               f   Z           �      �      �   �                  {   �             �   <   a   �   �   -   *   `   g   n   �   0           ?                   T          z           �       L      �       ,       3   x       �   �   �   �   d   �   j       I   �      [   �   �   9   �   2   7      8   W              U   ;      =   P           v      Jk     Qk  �  �  �  �  �    Xk         
   ����ok  2          �����k  0               �����k            �����k  <          ����5l  =               ���� 
  --delta[=OPTS]      Delta filter; valid OPTS (valid values; default):
                        dist=NUM   distance between bytes being subtracted
                                   from each other (1-256; 1) 
  --lzma1[=OPTS]      LZMA1 or LZMA2; OPTS is a comma-separated list of zero or
  --lzma2[=OPTS]      more of the following options (valid values; default):
                        preset=PRE reset options to a preset (0-9[e])
                        dict=NUM   dictionary size (4KiB - 1536MiB; 8MiB)
                        lc=NUM     number of literal context bits (0-4; 3)
                        lp=NUM     number of literal position bits (0-4; 0)
                        pb=NUM     number of position bits (0-4; 2)
                        mode=MODE  compression mode (fast, normal; normal)
                        nice=NUM   nice length of a match (2-273; 64)
                        mf=NAME    match finder (hc3, hc4, bt2, bt3, bt4; bt4)
                        depth=NUM  maximum search depth; 0=automatic (default) 
  --x86[=OPTS]        x86 BCJ filter (32-bit and 64-bit)
  --powerpc[=OPTS]    PowerPC BCJ filter (big endian only)
  --ia64[=OPTS]       IA-64 (Itanium) BCJ filter
  --arm[=OPTS]        ARM BCJ filter (little endian only)
  --armthumb[=OPTS]   ARM-Thumb BCJ filter (little endian only)
  --sparc[=OPTS]      SPARC BCJ filter
                      Valid OPTS for all BCJ filters:
                        start=NUM  start offset for conversions (default=0) 
 Basic file format and compression options:
 
 Custom filter chain for compression (alternative for using presets): 
 Operation modifiers:
 
 Other options:
 
With no FILE, or when FILE is -, read standard input.
       --block-list=SIZES
                      start a new .xz block after the given comma-separated
                      intervals of uncompressed data       --block-size=SIZE
                      start a new .xz block after every SIZE bytes of input;
                      use this to set the block size for threaded compression       --flush-timeout=TIMEOUT
                      when compressing, if more than TIMEOUT milliseconds has
                      passed since the previous flush and reading more input
                      would block, all pending data is flushed out       --ignore-check  don't verify the integrity check when decompressing       --info-memory   display the total amount of RAM and the currently active
                      memory usage limits, and exit       --memlimit-compress=LIMIT
      --memlimit-decompress=LIMIT
  -M, --memlimit=LIMIT
                      set memory usage limit for compression, decompression,
                      or both; LIMIT is in bytes, % of RAM, or 0 for defaults       --no-adjust     if compression settings exceed the memory usage limit,
                      give an error instead of adjusting the settings downwards       --no-sparse     do not create sparse files when decompressing
  -S, --suffix=.SUF   use the suffix `.SUF' on compressed files
      --files[=FILE]  read filenames to process from FILE; if FILE is
                      omitted, filenames are read from the standard input;
                      filenames must be terminated with the newline character
      --files0[=FILE] like --files but use the null character as terminator       --robot         use machine-parsable messages (useful for scripts)       --single-stream decompress only the first stream, and silently
                      ignore possible remaining input data       CheckVal %*s Header  Flags        CompSize    MemUsage  Filters   -0 ... -9           compression preset; default is 6; take compressor *and*
                      decompressor memory usage into account before using 7-9!   -F, --format=FMT    file format to encode or decode; possible values are
                      `auto' (default), `xz', `lzma', and `raw'
  -C, --check=CHECK   integrity check type: `none' (use with caution),
                      `crc32', `crc64' (default), or `sha256'   -Q, --no-warn       make warnings not affect the exit status   -T, --threads=NUM   use at most NUM threads; the default is 1; set to 0
                      to use as many threads as there are processor cores   -V, --version       display the version number and exit   -e, --extreme       try to improve compression ratio by using more CPU time;
                      does not affect decompressor memory requirements   -h, --help          display the short help (lists only the basic options)
  -H, --long-help     display this long help and exit   -h, --help          display this short help and exit
  -H, --long-help     display the long help (lists also the advanced options)   -k, --keep          keep (don't delete) input files
  -f, --force         force overwrite of output file and (de)compress links
  -c, --stdout        write to standard output and don't delete input files   -q, --quiet         suppress warnings; specify twice to suppress errors too
  -v, --verbose       be verbose; specify twice for even more verbose   -z, --compress      force compression
  -d, --decompress    force decompression
  -t, --test          test compressed file integrity
  -l, --list          list information about .xz files   Blocks:
    Stream     Block      CompOffset    UncompOffset       TotalSize      UncompSize  Ratio  Check   Blocks:             %s
   Check:              %s
   Compressed size:    %s
   Memory needed:      %s MiB
   Minimum XZ Utils version: %s
   Number of files:    %s
   Ratio:              %s
   Sizes in headers:   %s
   Stream padding:     %s
   Streams:
    Stream    Blocks      CompOffset    UncompOffset        CompSize      UncompSize  Ratio  Check      Padding   Streams:            %s
   Uncompressed size:  %s
  Operation mode:
 %s MiB of memory is required. The limit is %s. %s MiB of memory is required. The limiter is disabled. %s file
 %s files
 %s home page: <%s>
 %s:  %s: Cannot remove: %s %s: Cannot set the file group: %s %s: Cannot set the file owner: %s %s: Cannot set the file permissions: %s %s: Closing the file failed: %s %s: Error reading filenames: %s %s: Error seeking the file: %s %s: File already has `%s' suffix, skipping %s: File has setuid or setgid bit set, skipping %s: File has sticky bit set, skipping %s: File is empty %s: File seems to have been moved, not removing %s: Filename has an unknown suffix, skipping %s: Filter chain: %s
 %s: Input file has more than one hard link, skipping %s: Invalid argument to --block-list %s: Invalid filename suffix %s: Invalid multiplier suffix %s: Invalid option name %s: Invalid option value %s: Is a directory, skipping %s: Is a symbolic link, skipping %s: Not a regular file, skipping %s: Null character found when reading filenames; maybe you meant to use `--files0' instead of `--files'? %s: Options must be `name=value' pairs separated with commas %s: Read error: %s %s: Seeking failed when trying to create a sparse file: %s %s: Too many arguments to --block-list %s: Too small to be a valid .xz file %s: Unexpected end of file %s: Unexpected end of input when reading filenames %s: Unknown file format type %s: Unsupported integrity check type %s: Value is not a non-negative decimal integer %s: With --format=raw, --suffix=.SUF is required unless writing to stdout %s: Write error: %s %s: poll() failed: %s --list does not support reading from standard input --list works only on .xz files (--format=xz or --format=auto) 0 can only be used as the last element in --block-list Adjusted LZMA%c dictionary size from %s MiB to %s MiB to not exceed the memory usage limit of %s MiB Adjusted the number of threads from %s to %s to not exceed the memory usage limit of %s MiB Cannot establish signal handlers Cannot read data from standard input when reading filenames from standard input Compressed data cannot be read from a terminal Compressed data cannot be written to a terminal Compressed data is corrupt Compression and decompression with --robot are not supported yet. Compression support was disabled at build time Decompression support was disabled at build time Decompression will need %s MiB of memory. Disabled Empty filename, skipping Error creating a pipe: %s Error getting the file status flags from standard input: %s Error getting the file status flags from standard output: %s Error restoring the O_APPEND flag to standard output: %s Error restoring the status flags to standard input: %s Failed to enable the sandbox File format not recognized Internal error (bug) LZMA1 cannot be used with the .xz format Mandatory arguments to long options are mandatory for short options too.
 Maximum number of filters is four Memory usage limit for compression:     Memory usage limit for decompression:   Memory usage limit is too low for the given filter setup. Memory usage limit reached No No integrity check; not verifying file integrity None Only one file can be specified with `--files' or `--files0'. Report bugs to <%s> (in English or Finnish).
 Strms  Blocks   Compressed Uncompressed  Ratio  Check   Filename Switching to single-threaded mode due to --flush-timeout THIS IS A DEVELOPMENT VERSION NOT INTENDED FOR PRODUCTION USE. The .lzma format supports only the LZMA1 filter The environment variable %s contains too many arguments The exact options of the presets may vary between software versions. The filter chain is incompatible with --flush-timeout The sum of lc and lp must not exceed 4 Total amount of physical memory (RAM):  Totals: Try `%s --help' for more information. Unexpected end of input Unknown error Unknown-11 Unknown-12 Unknown-13 Unknown-14 Unknown-15 Unknown-2 Unknown-3 Unknown-5 Unknown-6 Unknown-7 Unknown-8 Unknown-9 Unsupported LZMA1/LZMA2 preset: %s Unsupported filter chain or filter options Unsupported options Unsupported type of integrity check; not verifying file integrity Usage: %s [OPTION]... [FILE]...
Compress or decompress FILEs in the .xz format.

 Using a preset in raw mode is discouraged. Valid suffixes are `KiB' (2^10), `MiB' (2^20), and `GiB' (2^30). Writing to standard output failed Yes Project-Id-Version: xz 5.2.6
Report-Msgid-Bugs-To: lasse.collin@tukaani.org
PO-Revision-Date: 2022-06-30 13:19+0200
Last-Translator: Remus-Gabriel Chelu <remusgabriel.chelu@disroot.org>
Language-Team: Romanian <translation-team-ro@lists.sourceforge.net>
Language: ro
MIME-Version: 1.0
Content-Type: text/plain; charset=UTF-8
Content-Transfer-Encoding: 8bit
X-Bugs: Report translation errors to the Language-Team address.
Plural-Forms: nplurals=3; plural=(n==1 ? 0 : n==0 || (n!=1 && n%100>=1 && n%100<=19) ? 1 : 2);
X-Generator: Poedit 2.3.1
 
  --delta[=OPȚIUNI]   OPȚIUNI valide (valori valide; valoarea implicită):
                        dist=NUM   distanța dintre octeți fiind
                                   scăzută unul de la celălalt (1-256; 1) 
  --lzma1[=OPȚIUNI]   LZMA1 sau LZMA2; OPȚIUNI este o listă separată prin
  --lzma2[=OPȚIUNI]   virgulă, de niciuna sau de mai multe dintre următoarele
                      opțiuni (între paranteze: valorile valide, și cele
                      predefinite):
                        preset=PRE resetează opțiunile la o presetare (0-9[e])
                        dict=NUM   dimensiunea dicționarului (4Ko - 1536Mo;
                                   8Mo)
                        lc=NUM     numărul de biți de context literal (0-4; 3)
                        lp=NUM     numărul de biți de poziție literală (0-4; 0)
                        pb=NUM     numărul de biți de poziție (0-4; 2)
                        mode=MOD   modul de comprimare (rapid, normal; normal)
                        nice=NUM   lungimea „drăguță” a unei potriviri (2-273;
                                   64)
                        mf=NUME    găsitor de potriviri (hc3, hc4, bt2, bt3,
                                   bt4; bt4)
                        depth=NUM  adâncimea maximă de căutare; 0=automată
                                   (valoarea implicită) 
  --x86[=OPȚIUNI]     filtru BCJ x86 (32-biți și 64-biți)
  --powerpc[=OPȚ]     filtru BCJ PowerPC (numai big endian)
  --ia64[=OPȚIUNI]    filtru BCJ IA-64 (Itanium)
  --arm[=OPȚIUNI]     filtru BCJ ARM  (numai little endian)
  --armthumb[=OPȚ]    filtru BCJ ARM-Thumb (numai little endian)
  --sparc[=OPȚIUNI]   filtru BCJ SPARC
                      OPȚIUNI valide pentru toate filtrele BCJ:
                        start=NUM poziția de pornire a conversiilor
                        (predefinit=0) 
 Formatul de bază al fișierului și opțiunile de comprimare:
 
 Lanț de filtre personalizat pentru comprimare (alternativă la utilizarea
 presetărilor): 
 Modificatori de operare:
 
 Alte opțiuni:
 
Fără FIȘIER, sau când FIȘIER este -, citește intrarea standard.
       --block-list=DM
                      începe un nou bloc .xz după intervalele date separate
                      prin virgulă, de date necomprimate       --block-size=DIM
                      începe un nou bloc .xz după fiecare DIM octeți de
                      intrare; utilizați acest lucru pentru a seta dimensiunea
                      blocului pentru compresia filară       --flush-timeout=TIMP_LIMITĂ
                      la comprimare, dacă au trecut mai mult de TIMP_LIMITĂ
                      milisecunde de la curățarea anterioară și citirea mai
                      multor intrări s-ar bloca, toate datele în așteptare
                      sunt eliminate       --ignore-check  nu se efectuează verificarea integrității la decomprimare       --info-memory   afișează cantitatea totală de memorie RAM și limitele de
                      utilizare a memoriei active în prezent, și iese       --memlimit-compress=LIMITĂ
      --memlimit-decompress=LIMITĂ
  -M, --memlimit=LIMITĂ
                      stabilește limita de utilizare a memoriei pentru
                      comprimare, decomprimare sau ambele; LIMITĂ este
                      exprimată în octeți, % din RAM sau 0 ca valori
                      predefinite       --no-adjust     dacă setările de comprimare depășesc limita de utilizare
                      a memoriei, dă o eroare în loc să reducă valorile setate       --no-sparse     nu creează fișiere împrăștiate când decomprimă
  -S, --suffix=.SUF   utilizează sufixul „.SUF” pentru fișierele comprimate
      --files[=FIȘIER]
                      citește numele fișierelor de procesat din FIȘIER;
                      dacă FIȘIER este omis, numele de fișiere sunt citite
                      de la intrarea standard; numele de fișiere trebuie să
                      fie terminate cu caracterul de linie\ nouăn
      --files0[=FIȘIER]
                      ca --files, dar folosește caracterul nul(null) ca
                      terminator       --robot         utilizează mesaje analizabile de mașină (utile
                      pentru scripturi)       --single-stream decomprimă doar primul flux și ignoră în tăcere
                      posibilele date de intrare rămase       Val.Verif. %*sAntet  Indicat.    Dim.Comp.      UzMem.  Filtre   -0 ... -9           presetare comprimare; valoarea predefinită este 6;
                      luați în considerare memoria utilizată de compresor
                      *și* decompresor, înainte de a utiliza valorile 7-9!   -F, --format=FORMAT
                      formatul de fișier pentru codificare sau decodare;
                      valorile posibile sunt „auto” (predefinit), „xz”,
                      „lzma” și „raw”
  -C, --check=VERIFICARE
                      tip de verificare a integrității: „niciuna” (utilizați
                      cu precauție), „crc32”, „crc64” (predefinit) sau
                      „sha256”   -Q, --no-warn       face ca avertismentele să nu afecteze starea de ieșire   -T, --threads=NUMĂR
                      utilizați cel mult NUMĂR fire de execuție; predefinit
                      este 1; setați la 0 pentru a utiliza atâtea fire de
                      execuție câte nuclee de procesor există   -V, --version       afișează numărul versiunii, și iese   -e, --extreme       încearcă să îmbunătățească raportul de comprimare
                      folosind mai mult timp CPU-ul; nu afectează cerințele
                      de memorie ale decompresorului   -h, --help          afișează ajutorul scurt (enumeră doar opțiunile de bază)
  -H, --long-help     afișează acest ajutor lung(detaliat) și iese   -h, --help          afișează acest scurt mesaj de ajutor, și iese
  -H, --long-help     afișează mesajul detaliat de ajutor (afișează și opțiunile
                      avansate)   -k, --keep          păstrează (nu șterge) fișierele de intrare
  -f, --force         forțează suprascrierea fișierului de ieșire și
                      (de)comprimă legăturile
  -c, --stdout        scrie la ieșirea standard și nu șterge fișierele de
                      intrare   -q, --quiet         suprimă avertismentele; specificați de două ori pentru
                      a suprima și erorile
  -v, --verbose       descriere detaliată; specificați de două ori pentru și
                      mai multe detalii   -z, --compress      forțează comprimarea
  -d, --decompress    forțează decomprimarea
  -t, --test          testează integritatea fișierului comprimat
  -l, --list          listează informații despre fișierele .xz   Blocuri:
      Flux      Bloc    DecalajComp.  DecalajDecomp.      Dim.Totală     Dim.Decomp. Raport Verif.   Blocuri:                   %s
   Verificare:                %s
   Dimensiune comprimată:     %s
   Memorie necesară:          %s MiO
   Versiunea minimă XZ Utils: %s
   Numărul de fișiere:        %s
   Raport:                    %s
   Dim. în antete:            %s
   Completare flux:           %s
   Fluxuri:
      Flux   Blocuri    DecalajComp.  DecalajDecomp.       Dim.Comp.     Dim.Decomp. Raport Verif.   Compl.Flux   Fluxuri:                   %s
   Dimensiune decomprimată:   %s
  Mod de operare:
 Se necesită %s MiO de memorie. Limita este de %s. Se necesită %s MiO de memorie. Limitatorul este dezactivat. %s fișier
 %s fișiere
 %s de fișiere
 %s pagina principală: <%s>
 %s:  %s: Nu se poate elimina: %s %s: Nu se poate configura grupul proprietar al fișierului: %s %s: Nu se poate configura proprietarul fișierului: %s %s: Nu se pot configura permisiunile fișierului: %s %s: Închiderea fișierului a eșuat: %s %s: Eroare la citirea numelor de fișiere: %s %s: Eroare la căutarea fișierului: %s %s: Fișierul are deja sufixul „%s”, se omite %s: Fișierul are activați biții «setuid» sau «setgid», se omite %s: Fișierul are activat bitul lipicios(sticky), se omite %s: Fișierul este gol %s: Fișierul pare să fi fost mutat, nu eliminat %s: Numele fișierului are un sufix necunoscut, care se omite %s: Lanț de filtre: %s
 %s: Fișierul de intrare are mai mult de o legătură dură, se omite %s: Argument nevalid pentru opțiunea „--block-list” %s: Sufixul numelui de fișier nu este valid %s: Sufix multiplicator nevalid %s: Nume de opțiune nevalid %s: Valoare nevalidă a opțiunii %s: Este un director, se omite %s: Este o legătură simbolică, se omite %s: Nu este un fișier obișnuit, se omite %s: Caracter nul(null) găsit la citirea numelor de fișiere; poate ați vrut să utilizați „--files0” în loc de „--files”? %s: Opțiunile trebuie să fie perechi „nume=valoare” separate prin virgule %s: Eroare de citire: %s %s: Căutarea a eșuat când se încerca crearea unui fișier împrăștiat, disipat(sparse): %s %s: Prea multe argumente pentru opțiunea „--block-list” %s: Prea mic pentru a fi un fișier .xz valid %s: Sfârșit neașteptat al fișierului %s: Sfârșit neașteptat al intrării la citirea numelor de fișiere %s: Tip de format de fișier necunoscut %s: Tip de verificare a integrității neacceptat %s: Valoarea nu este un număr întreg zecimal nenegativ %s: Cu --format=raw, este necesar --suffix=.SUF, cu excepția cazului în care se scrie la ieșirea standard(stdout) %s: Eroare de scriere: %s %s: poll() a eșuat: %s --list nu acceptă citirea de la intrarea standard --list funcționează numai pe fișierele .xz (--format=xz sau --format=auto) 0 poate fi folosit doar ca ultimul element din opțiunea „--block-list” S-a ajustat dimensiunea dicționarului LZMA%c de la %s MiO la %s MiO pentru a nu se depăși limita de utilizare a memoriei de %s MiO Numărul de fire de execuție a fost ajustat de la %s la %s pentru a nu se depăși limita de utilizare a memoriei de %s MiO Nu se pot stabili operatorii de semnal Nu se pot citi date de la intrarea standard atunci când se citesc numele de fișiere de la intrarea standard Datele comprimate nu pot fi citite de pe un terminal Datele comprimate nu pot fi scrise pe un terminal Datele comprimate sunt corupte Comprimarea și decomprimarea cu --robot nu sunt încă acceptate. Suportul de comprimare a fost dezactivat în timpul construirii Suportul de decomprimare a fost dezactivat în timpul construirii Decompresia va avea nevoie de %s MiO de memorie. Dezactivat Numele fișierului este gol, se omite Eroare la crearea unui racord(pipe): %s Eroare la obținerea indicatorilor de stare a fișierului de la intrarea standard: %s Eroare la obținerea indicatorilor de stare a fișierului de la ieșirea standard: %s Eroare la restabilirea indicatorului O_APPEND la ieșirea standard: %s Eroare la restabilirea indicatorilor de stare la intrarea standard: %s Nu s-a reușit activarea cutiei de probă(sandbox) Formatul fișierului nu este recunoscut Eroare internă (bug) LZMA1 nu poate fi utilizat cu formatul .xz Argumentele obligatorii pentru opțiunile lungi sunt obligatorii și
pentru opțiunile scurte.
 Numărul maxim de filtre este patru Limita de utilizare a memoriei pentru comprimare:   Limita de utilizare a memoriei pentru decomprimare: Limita de utilizare a memoriei este prea mică pentru configurarea dată filtrului. Limita de utilizare a memoriei a fost atinsă Nu Fără verificare a integrității; nu se verifică integritatea fișierului Nici una Numai un fișier poate fi specificat cu „--files” sau „--files0”. Raportați erorile la <%s> (în engleză sau finlandeză).
Raportați erorile de traducere la <translation-team-ro@lists.sourceforge.net>
 Fluxuri Blocuri Comprimare Decomprimare Raport  Verificare Nume fișier Se trece la modul cu un singur fir datorită opțiunii „--flush-timeout” ACEASTA ESTE O VERSIUNE DE DEZVOLTARE, NEDESTINATĂ UTILIZĂRII ÎN PRODUCȚIE. Formatul .lzma acceptă numai filtrul LZMA1 Variabila de mediu „%s” conține prea multe argumente Opțiunile exacte ale presetărilor pot varia între versiunile de software. Lanțul de filtre este incompatibil cu opțiunea „--flush-timeout” Suma de lc și lp nu trebuie să depășească 4 Cantitatea totală de memorie fizică (RAM):          Totaluri: Încercați «%s --help» pentru mai multe informații. Sfârșit neașteptat al intrării Eroare necunoscută Necunos-11 Necunos-12 Necunos-13 Necunos-14 Necunos-15 Necunos-2 Necunos-3 Necunos-5 Necunos-6 Necunos-7 Necunos-8 Necunos-9 Presetare LZMA1/LZMA2 neacceptată: %s Lanț de filtre sau opțiuni de filtrare neacceptate Opțiuni neacceptate Tip neacceptat de verificare a integrității; nu se verifică integritatea fișierului Utilizare: %s [OPȚIUNE]... [FIȘIER]...
Comprimă sau decomprimă FIȘIER(e) în formatul .xz.

 Utilizarea unei presetări în modul brut este descurajată. Sufixele valide sunt „KiO” (2^10), „MiO” (2^20) și „GiO” (2^30). Scrierea la ieșirea standard a eșuat Da PRIu32 PRIu64 Using up to % threads. The selected match finder requires at least nice=% Value of the option `%s' must be in the range [%, %] Se utilizează până la % fire. Găsitorul de potriviri selectat necesită cel puțin nice=% Valoarea opțiunii „%s” trebuie să fie în intervalul [%, %] 