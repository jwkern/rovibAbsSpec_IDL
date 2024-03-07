;Working to clean up the code and delete superfluous lines of test code and other
;miscellaneous plots
;
;Update likely March 6th-8th
;
;Create synthetic spectrum of optically thick LTE gas 
;assuming gaussian lineshape. 
; i.e. -A0*exp(-(freq-A1)^2/(A2^2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;; DEFINITIONS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FOR iii=0,1 DO BEGIN
	X1213   = 50.0000
	X1218   = 500.
	inst_res=3.50
	f=fbig
	ratio=rbig

	;Define constants
	v	=double(1)		;velocity resolution in km/s
	c	=double(2.997924562E5)	;speed of light in km/s
	hc	=1.986484E-16		;erg*cm
	k	=1.380622E-16		;erg/K
	kwn     =0.6950356              ;cm-1/K
	B	=1.92			; cm-1
	me      =9.1094d-28
	elec_ch =4.8032d-10


	;;;;;;;;;;;;;;;; DOPPLER SHIFTS of COMPONENTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	ds12_a=-21.8
	ds12_b=-4.00
	ds12_bb=-1.8
	ds12_cc=4.0
	ds12_c=5.5  
	ds12_d=12.5

	ds1221=10.2

	ds13_b=-4.0
	ds13_bb=-1.8
	ds13_cc=4.25
	ds13_c=5.50
	ds13_d=12.5

	ds18_bb=5.8 
	ds18_c=ds12_c
	ds18_cc=ds12_cc
	ds18_d=ds12_d

	ds1221cc=ds12_cc
	ds1221c=ds12_c
	ds1221d=ds12_d


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;Define Freq Scale
	freq=[2190.d]
	i=1
	REPEAT BEGIN
	   freq=[freq,freq[i-1]-freq[i-1]*v/c]
	   i=i+1
	ENDREP UNTIL MIN(freq) LE 1990.0

	;freq=DINDGEN(30720,start=1900,increment=0.013020833)



	;;;;;;;;;;;;;;;;;;;;;;; CO MOLECULAR PROPERTIES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;READCOL,'/home/jwkern/Research/Exopl/Processed_Data/IRTF/ishell/2019/02feb/Models/CO_hitran_v10.txt',line_freq12,EinA12,Eup12,Elow12,vup12,vlow12,Jup12,Jlow12,gup12,glow12
	;READCOL,'/home/jwkern/Research/Exopl/Processed_Data/IRTF/ishell/2019/02feb/Models/CO_hitran_v21.txt',line_freq1221,EinA1221,Eup1221,Elow1221,vup1221,vlow1221,Jup1221,Jlow1221,gup1221,glow1221
	;READCOL,'/home/jwkern/Research/Exopl/Processed_Data/IRTF/ishell/2019/02feb/Models/13CO_hitran_v10.txt',line_freq13,EinA13,Eup13,Elow13,vup13,vlow13,Jup13,Jlow13,gup13,glow13
	;READCOL,'/home/jwkern/Research/Exopl/Processed_Data/IRTF/ishell/2019/02feb/Models/18CO_hitran_v10.txt',line_freq18,EinA18,Eup18,Elow18,vup18,vlow18,Jup18,Jlow18,gup18,glow18

	osc_str12=me*c*1e5*gup12*EinA12/(8.*!pi^2*elec_ch^2*line_freq12^2*glow12)
	osc_str1221=me*c*1e5*gup1221*EinA1221/(8.*!pi^2*elec_ch^2*line_freq1221^2*glow1221)
	osc_str13=me*c*1e5*gup13*EinA13/(8.*!pi^2*elec_ch^2*line_freq13^2*glow13)
	osc_str18=me*c*1e5*gup18*EinA18/(8.*!pi^2*elec_ch^2*line_freq18^2*glow18)



	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;;;;;;;;;;;;;;;;;;;;; CALCULATE MODEL SPECTRUM ;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	size=3 

	;window,3
	;window,2
	;window,1
	;window,0

	;FOR j=0,size-1 DO BEGIN         
	veil_a=0.20
	veil_cc=0.8;1.05
	veil_cc13=0.08
	veil_c =0.60
	veil_c13=0.6
	veil_d=0.3
	veil_d13=0.3

	veil_v21_d=0.3
        veil_v21_c=0.6
        veil_v21_cc=0.8


	;;;;;;;;;;;;;;;;;;;;;;; 12CO v=1-0 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
	trot12a    =991;RANDOMU(seed)*600.0 + 900 ;900.11    
	N12a       =8.0559001e17;10.^(RANDOMU(seed)*0.3 + ALOG10(2.3939000e+17) - 0.15)
	int_wth12a =2.1;9.595081;RANDOMU(seed)*2.0 + 11.130014 - 1.0

	trot12b    =126.0700;RANDOMU(seed)*100.0 + 150.46490 - 50.0
	N12b       =5.0003000e+17;10.^(RANDOMU(seed)*0.3 + ALOG10(1.2335000e+17) - 0.15)
	int_wth12b =2.1522020;RANDOMU(seed)*2.0 + 3.3340890 - 1.0

	;trot12bb   =RANDOMU(seed)*150.0 + 482.20639 - 75.0
	;N12bb      =10.^(RANDOMU(seed)*0.3 + ALOG10(1.4280000e+17) - 0.15)
	;int_wth12bb=RANDOMU(seed)*2.0 + 5.3606339 - 1.0

	
	trot12cc   =276.6;RANDOMU(seed)*400.0 + 200.6 - 150.0
	N12cc      =5.00e+18;10.^(RANDOMU(seed)*0.3 + ALOG10(2.6407000e+17) - 0.15)
	int_wth12cc=2.10;RANDOMU(seed)*2.0 + 8.8413506 - 1.0


        trot12c    =500;DOMU(seed)*1000.0 + 700.0 - 500.0;1174.3561 - 100.0
        N12c       =7.00e18;10.^(RANDOMU(seed)*1.0 + ALOG10(6.376500e+17) - 0.15) 
        int_wth12c =2.15;RANDOMU(seed)*4.0 + 5.2684239 - 4.0



	trot12d    =830.0;ANDOMU(seed)*200.0 + 1494.1101 - 100.00 
	N12d       =9.30e+18;10.^(RANDOMU(seed)*0.3 + ALOG10(3.2577001e+17) - 0.15)
	int_wth12d =2.21;ANDOMU(seed)*2.0 + 5.8233671 - 1.0





	;;;;;;;;;;;;;;;;;;;;;;; 12CO v=2-1 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
	c1=10
	
	trot1221cc    =trot12cc
	N1221cc       =N12cc/11
	int_wth1221cc =int_wth12cc

	trot1221c    =trot12c
        N1221c       =N12c/11 
        int_wth1221c =int_wth12c

	trot1221d    =trot12d
        N1221d       =N12d/11
        int_wth1221d =int_wth12d




	;;;;;;;;;;;;;;;;;;;;;;; 13CO v=1-0 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	trot13b    =108;trot12b
	X1213b	   =50;(RANDOMU(seed)*10 + 45)
	;print,X1213b
	N13b       =N12b/X1213b;1.1599999e+16
	int_wth13b =int_wth12b

	;trot13bb    =110.00
	;N13bb       =3.2755000e+16
	;int_wth13bb =int_wth12bb

	trot13cc    =201;trot12cc
	X1213cc     =30.;(RANDOMU(seed)*30 + 5)
	;print,X1213cc
	N13cc       =N12cc/X1213cc;3.5199001e+16
	int_wth13cc =int_wth12cc-0.20


	trot13c    =407;trot12c
        X1213c     =85;(RANDOMU(seed)*10 + 45)
        ;print,X1213c
        N13c       =N12c/X1213c;35;7.1199001e+17
        int_wth13c =int_wth12c-0.25




	trot13d    =702;12d
	X1213d     =87;(RANDOMU(seed)*10 + 45)
	;print,X1213d
	N13d       =N12d/X1213d;1.6411000e+17
	int_wth13d =int_wth13d-0.35



	;;;;;;;;;;;;;;;;;;;;;;; C18O v=1-0 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 


	;trot18bb     =425.0 ;RANDOMU(seed)*30. + 260.67 - 15.
	;N18bb        =0.8624e+15;10.^(RANDOMU(seed)*0.1 + 16.0) + 2.8624000e+16 - 1.5848932e+16
	;int_wth18bb  =2.6862190;RANDOMU(seed)*0.22 + 2.5 

	trot18c       =421;trot12c;325.0 ;RANDOMU(seed)*30. + 260.67 - 15.
	N18c          =N12c/450.0;3.1024e+16;10.^(RANDOMU(seed)*0.1 + 16.0) + 2.8624000e+16 - 1.5848932e+16
	int_wth18c    =int_wth12c-0.25;RANDOMU(seed)*0.22 + 2.5 

	trot18cc    =203;trot12cc;130.0;RANDOMU(seed)*30. + 260.67 - 15.
	N18cc       =N12cc/405.0;4.2624e+15;10.^(RANDOMU(seed)*0.1 + 16.0) + 2.8624000e+16 - 1.5848932e+16
	int_wth18cc =int_wth12cc-0.25;2190;RANDOMU(seed)*0.22 + 2.5

	trot18d     =trot12d;850.0 ;RANDOMU(seed)*30. + 260.67 - 15.
	N18d        =N12d/550.0;3.0624e+16;10.^(RANDOMU(seed)*0.1 + 16.0) + 2.8624000e+16 - 1.5848932e+16
	int_wth18d  =int_wth12d;6.00;RANDOMU(seed)*0.22 + 2.5


	;goto,step_COG

	;determine values of the partition function, state densities, and optical
	;depths to setup variables needed for radiative transfer eq. 
	index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot12a)),1))
	Q12a=Q_T12(index)
	;Q12a=part_CO(trot12a)
	NJ12a=N12a*glow12*exp(-Elow12/(kwn*trot12a))/Q12a(0)
	tau_0_12a=double(NJ12a*osc_str12/(66.8*int_wth12a*1e5*line_freq12))
	print,MAX(tau_0_12a)

	index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot12b)),1))
	Q12b=Q_T12(index)
	;Q12b=part_CO(trot12b)
	NJ12b=N12b*glow12*exp(-Elow12/(kwn*trot12b))/Q12b(0)
	;print,NJ12b(0)
	tau_0_12b=double(NJ12b*osc_str12/(66.8*int_wth12b*1e5*line_freq12))
	print,MAX(tau_0_12b)

	;index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot12bb)),1))
	;Q12bb=Q_T12(index)
	;Q12b=part_CO(trot12b)
	;NJ12bb=N12bb*gup12*exp(-Eup12/(kwn*trot12bb))/Q12bb(0)
	;print,NJ12b(0)
	;tau_0_12bb=double(NJ12bb*osc_str12/(66.8*int_wth12bb*1e5*line_freq12))


	index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot12c)),1))
	Q12c=Q_T12(index)
	;Q12c=part_CO(trot12c)
	NJ12c=N12c*glow12*exp(-Elow12/(kwn*trot12c))/Q12c(0)
	tau_0_12c=double(NJ12c*osc_str12/(66.8*int_wth12c*1e5*line_freq12))
	print,MAX(tau_0_12c)



	index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot12cc)),1))
	Q12cc=Q_T12(index)
	;Q12c=part_CO(trot12c)
	NJ12cc=N12cc*glow12*exp(-Elow12/(kwn*trot12cc))/Q12cc(0)
	tau_0_12cc=double(NJ12cc*osc_str12/(66.8*int_wth12cc*1e5*line_freq12))
	print,MAX(tau_0_12cc)


	index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot12d)),1))
	Q12d=Q_T12(index)
	;Q12d=part_CO(trot12d)
	NJ12d=N12d*glow12*exp(-Elow12/(kwn*trot12d))/Q12d(0)
	tau_0_12d=double(NJ12d*osc_str12/(66.8*int_wth12d*1e5*line_freq12))
	print,MAX(tau_0_12d)


	index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot1221cc)),1))
	Q1221cc=Q_T12(index)
	;Q1221=part_CO(trot1221)
	NJ1221cc=N1221cc*glow1221*exp(-Elow1221/(kwn*trot1221cc))/Q1221cc(0)
	tau_0_1221cc=double(NJ1221cc*osc_str1221/(66.8*int_wth1221cc*1e5*line_freq1221))
	print,MAX(tau_0_1221cc)



	index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot1221c)),1))
        Q1221c=Q_T12(index)
        ;Q1221=part_CO(trot1221)
        NJ1221c=N1221c*glow1221*exp(-Elow1221/(kwn*trot1221c))/Q1221c(0)
        tau_0_1221c=double(NJ1221c*osc_str1221/(66.8*int_wth1221c*1e5*line_freq1221))
        print,MAX(tau_0_1221c)


	index=WHERE(T_Q12 EQ STRTRIM(STRING(FIX(trot1221d)),1))
        Q1221d=Q_T12(index)
        ;Q1221=part_CO(trot1221)
        NJ1221d=N1221d*glow1221*exp(-Elow1221/(kwn*trot1221d))/Q1221d(0)
        tau_0_1221d=double(NJ1221d*osc_str1221/(66.8*int_wth1221d*1e5*line_freq1221))
        print,MAX(tau_0_1221d)




	index=WHERE(T_Q13 EQ STRTRIM(STRING(FIX(trot13b)),1))
	Q13b=Q_T13(index)
	;Q13a=part_13CO(trot13b)
	NJ13b=N13b*glow13*exp(-Elow13/(kwn*trot13b))/Q13b(0)
	tau_0_13b=double(NJ13b*osc_str13/(66.8*int_wth13b*1e5*line_freq13))
	print,MAX(tau_0_13b)




	;index=WHERE(T_Q13 EQ STRTRIM(STRING(FIX(trot13bb)),1))
	;Q13bb=Q_T13(index)
	;Q13b=part_13CO(trot13bb)
	;NJ13bb=N13bb*gup13*exp(-Eup13/(kwn*trot13bb))/Q13bb(0)
	;tau_0_13bb=double(NJ13bb*osc_str13/(66.8*int_wth13bb*1e5*line_freq13))

	index=WHERE(T_Q13 EQ STRTRIM(STRING(FIX(trot13c)),1))
	Q13c=Q_T13(index)
	;Q13c=part_13CO(trot13c)
	NJ13c=N13c*glow13*exp(-Elow13/(kwn*trot13c))/Q13c(0)
	tau_0_13c=double(NJ13c*osc_str13/(66.8*int_wth13c*1e5*line_freq13))
	print,MAX(tau_0_13c)




	index=WHERE(T_Q13 EQ STRTRIM(STRING(FIX(trot13cc)),1))
	Q13cc=Q_T13(index)
	;Q13c=part_13CO(trot13c)
	NJ13cc=N13cc*glow13*exp(-Elow13/(kwn*trot13cc))/Q13cc(0)
	tau_0_13cc=double(NJ13cc*osc_str13/(66.8*int_wth13cc*1e5*line_freq13))
	print,MAX(tau_0_13cc)


	index=WHERE(T_Q13 EQ STRTRIM(STRING(FIX(trot13d)),1))
	Q13d=Q_T13(index)
	;Q13d=part_13CO(trot13d)
	NJ13d=N13d*glow13*exp(-Elow13/(kwn*trot13d))/Q13d(0)
	tau_0_13d=double(NJ13d*osc_str13/(66.8*int_wth13d*1e5*line_freq13))
	print,MAX(tau_0_13d)



	;index=WHERE(T_Q18 EQ STRTRIM(STRING(FIX(trot18bb)),1))
	;Q18bb=Q_T18(index)
	;Q18=part_C18O(trot18bb)
	;NJ18bb=N18bb*gup18*exp(-Eup18/(kwn*trot18bb))/Q18bb(0)
	;tau_0_18bb=double(NJ18bb*osc_str18/(66.8*int_wth18bb*1e5*line_freq18))


	index=WHERE(T_Q18 EQ STRTRIM(STRING(FIX(trot18c)),1))
	Q18c=Q_T18(index)
	;Q18=part_C18O(trot18bb)
	NJ18c=N18c*glow18*exp(-Elow18/(kwn*trot18c))/Q18c(0)
	tau_0_18c=double(NJ18c*osc_str18/(66.8*int_wth18c*1e5*line_freq18))
	print,MAX(tau_0_18c)

	index=WHERE(T_Q18 EQ STRTRIM(STRING(FIX(trot18cc)),1))
	Q18cc=Q_T18(index)
	;Q18=part_C18O(trot18bb)
	NJ18cc=N18cc*glow18*exp(-Elow18/(kwn*trot18cc))/Q18cc(0)
	tau_0_18cc=double(NJ18cc*osc_str18/(66.8*int_wth18cc*1e5*line_freq18))
	print,MAX(tau_0_18cc)

	index=WHERE(T_Q18 EQ STRTRIM(STRING(FIX(trot18d)),1))
	Q18d=Q_T18(index)
	;Q18=part_C18O(trot18bb)
	NJ18d=N18d*glow18*exp(-Elow18/(kwn*trot18d))/Q18d(0)
	tau_0_18d=double(NJ18d*osc_str18/(66.8*int_wth18d*1e5*line_freq18))
	print,MAX(tau_0_18d)

	;Calculate tau(nu) for system components
	tau_diskabs_a=FLTARR(N_ELEMENTS((freq)))
	tau_cloud=tau_diskabs_a
	tau_cloud_cc=tau_diskabs_a
	tau_diskabs_cc13=tau_diskabs_a
	;tau_diskem12=tau_diskabs_a
	;tau_diskem13=tau_diskabs_a
	;tau_diskem18=tau_diskabs_a
	tau_diskabs_c=tau_diskabs_a
	tau_diskabs_c13=tau_diskabs_a
	tau_diskabs_d=tau_diskabs_a
	tau_diskabs_d13=tau_diskabs_a
	tau_diskabs_v21_d=tau_diskabs_a
        tau_diskabs_v21_c=tau_diskabs_a
        tau_diskabs_v21_cc=tau_diskabs_a




	;;;;;;;;;;;; DISK ABSORPTION COMPONENTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	FOR i=0,N_ELEMENTS(line_freq12)-1  DO tau_diskabs_a=tau_diskabs_a $
	   + INTERPOL(tau_0_12a(i)*exp(-(c*(double(1.)-double(freq/line_freq12(i)))/int_wth12a)^2),freq,freq+ds12_a*freq/c)

	FOR i=0,N_ELEMENTS(line_freq12)-1  DO tau_diskabs_c=tau_diskabs_c $
	   + INTERPOL(tau_0_12c(i)*exp(-(c*(double(1.)-double(freq/line_freq12(i)))/int_wth12c)^2),freq,freq+ds12_c*freq/c) 

	FOR i=0,N_ELEMENTS(line_freq12)-1  DO tau_diskabs_d=tau_diskabs_d $
	   + INTERPOL(tau_0_12d(i)*exp(-(c*(double(1.)-double(freq/line_freq12(i)))/int_wth12d)^2),freq,freq+ds12_d*freq/c)
	



	FOR i=0,N_ELEMENTS(line_freq1221)-1  DO tau_diskabs_v21_d=tau_diskabs_v21_d $
	   + INTERPOL(tau_0_1221d(i)*exp(-(c*(double(1.)-double(freq/line_freq1221(i)))/int_wth1221d)^2),freq,freq+ds1221d*freq/c)

        FOR i=0,N_ELEMENTS(line_freq1221)-1  DO tau_diskabs_v21_c=tau_diskabs_v21_c $
           + INTERPOL(tau_0_1221c(i)*exp(-(c*(double(1.)-double(freq/line_freq1221(i)))/int_wth1221c)^2),freq,freq+ds1221c*freq/c)

           FOR i=0,N_ELEMENTS(line_freq1221)-1  DO tau_diskabs_v21_cc=tau_diskabs_v21_cc $
           + INTERPOL(tau_0_1221cc(i)*exp(-(c*(double(1.)-double(freq/line_freq1221(i)))/int_wth1221cc)^2),freq,freq+ds1221cc*freq/c)



	FOR i=0,N_ELEMENTS(line_freq13)-1  DO tau_diskabs_c13=tau_diskabs_c13 $
	   + INTERPOL(tau_0_13c(i)*exp(-(c*(double(1.)-double(freq/line_freq13(i)))/int_wth13c)^2),freq,freq+ds13_c*freq/c) 

	FOR i=0,N_ELEMENTS(line_freq13)-1  DO tau_diskabs_d13=tau_diskabs_d13 $
	   + INTERPOL(tau_0_13d(i)*exp(-(c*(double(1.)-double(freq/line_freq13(i)))/int_wth13d)^2),freq,freq+ds13_d*freq/c)





	FOR i=0,N_ELEMENTS(line_freq18)-1  DO tau_diskabs_c13=tau_diskabs_c13 $
	   + INTERPOL(tau_0_18c(i)*exp(-(c*(double(1.)-double(freq/line_freq18(i)))/int_wth18c)^2),freq,freq+ds18_c*freq/c) 


	FOR i=0,N_ELEMENTS(line_freq18)-1  DO tau_diskabs_d13=tau_diskabs_d13 $ 
	   + INTERPOL(tau_0_18d(i)*exp(-(c*(double(1.)-double(freq/line_freq18(i)))/int_wth18d)^2),freq,freq+ds18_d*freq/c) 





	;;;;;;;;;; ENVELOPE ABSORPTION COMPONENTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	FOR i=0,N_ELEMENTS(line_freq12)-1 DO tau_cloud=tau_cloud + INTERPOL(tau_0_12b(i)*exp(-(c*(double(1.)-double(freq/line_freq12(i)))/int_wth12b)^2),freq,freq+ds12_b*freq/c)  

	FOR i=0,N_ELEMENTS(line_freq13)-1 DO tau_cloud=tau_cloud + INTERPOL(tau_0_13b(i)*exp(-(c*(double(1.)-double(freq/line_freq13(i)))/int_wth13b)^2),freq,freq+ds13_b*freq/c) 



	FOR i=0,N_ELEMENTS(line_freq18)-1 DO tau_diskabs_cc13=tau_diskabs_cc13 + INTERPOL(tau_0_18cc(i)*exp(-(c*(double(1.)-double(freq/line_freq18(i)))/int_wth18cc)^2),freq,freq+ds18_cc*freq/c)


	FOR i=0,N_ELEMENTS(line_freq12)-1 DO tau_cloud_cc=tau_cloud_cc + INTERPOL(tau_0_12cc(i)*exp(-(c*(double(1.)-double(freq/line_freq12(i)))/int_wth12cc)^2),freq,freq+ds12_cc*freq/c)

	FOR i=0,N_ELEMENTS(line_freq13)-1 DO tau_diskabs_cc13=tau_diskabs_cc13 + INTERPOL(tau_0_13cc(i)*exp(-(c*(double(1.)-double(freq/line_freq13(i)))/int_wth13cc)^2),freq,freq+ds13_cc*freq/c)




	;;;;;;;;;; DISK EMISSION COMPONENTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;FOR i=0,N_ELEMENTS(line_freq12)-1 DO tau_diskem12=tau_diskem12 + INTERPOL(tau_0_12bb(i)*exp(-(c*(double(1.)-double(freq/line_freq12(i)))/int_wth12bb)^2),freq,freq+ds12_bb*freq/c) 

	;FOR i=0,N_ELEMENTS(line_freq13)-1 DO tau_diskem13=tau_diskem13 + INTERPOL(tau_0_13bb(i)*exp(-(c*(double(1.)-double(freq/line_freq13(i)))/int_wth13bb)^2),freq,freq+ds13_bb*freq/c) 

	;FOR i=0,N_ELEMENTS(line_freq18)-1 DO tau_diskem18=tau_diskem18 + INTERPOL(tau_0_18bb(i)*exp(-(c*(double(1.)-double(freq/line_freq18(i)))/int_wth18bb)^2),freq,freq+ds18_bb*freq/c)

	;B_T12=2*hc*c*1e5*(freq^3)/(exp((hc*freq)/(k*trot12bb))-1)
	;B_T13=2*hc*c*1e5*(freq^3)/(exp((hc*freq)/(k*trot13bb))-1)
	;B_T18=2*hc*c*1e5*(freq^3)/(exp((hc*freq)/(k*trot18bb))-1)




	spectrum_disk_c=(((exp(-tau_diskabs_c)+veil_c)/(veil_c+1)))
	spectrum_disk_d=(((exp(-tau_diskabs_d)+veil_d)/(veil_d+1)))
	spectrum_wind=(((exp(-tau_diskabs_a)+veil_a)/(veil_a+1)))
	spectrum_env=exp(-tau_cloud)
	spectrum_disk_cc=(((exp(-tau_cloud_cc)+veil_cc)/(veil_cc+1)))
	spectrum_disk_cc13=(((exp(-tau_diskabs_cc13)+veil_cc13)/(veil_cc13+1)))
	spectrum_disk_c13=(((exp(-tau_diskabs_c13)+veil_c13)/(veil_c13+1)))
        spectrum_disk_d13=(((exp(-tau_diskabs_d13)+veil_d13)/(veil_d13+1)))
	spectrum_disk_v21_cc=(((exp(-tau_diskabs_v21_cc)+veil_v21_cc)/(veil_v21_cc+1)))
        spectrum_disk_v21_c=(((exp(-tau_diskabs_v21_c)+veil_v21_c)/(veil_v21_c+1)))
        spectrum_disk_v21_d=(((exp(-tau_diskabs_v21_d)+veil_v21_d)/(veil_v21_d+1)))








	disk_res=6.5;radiative transfer equation for given tau(nu)'s
	;currently setup with a disk emitting through a molecular cloud
	vfreq=(freq(9500)-freq)*c/freq(9500)
	a3=(disk_res)/1.665
	line_prof1=-exp(-(vfreq/a3)^2)/(SQRT(!pi)*a3)
	line_prof1=SHIFT(line_prof1,-9500)
	conv_spec=REAL_PART(FFT(FFT(spectrum_disk_c)*FFT(line_prof1),1))
	spectrum_disk_c=conv_spec*total(spectrum_disk_c)/total(conv_spec) ; <--SPECTRUM of DISK!!!
	conv_spec=REAL_PART(FFT(FFT(spectrum_disk_c13)*FFT(line_prof1),1))
        spectrum_disk_c13=conv_spec*total(spectrum_disk_c13)/total(conv_spec) ; <--SPECTRUM of DISK!!!
	conv_spec=REAL_PART(FFT(FFT(spectrum_disk_v21_c)*FFT(line_prof1),1))
        spectrum_disk_v21_c=conv_spec*total(spectrum_disk_v21_c)/total(conv_spec) ; <--SPECTRUM of DISK!!!




	disk2_res=12.0;radiative transfer equation for given tau(nu)'s
	;currently setup with a disk emitting through a molecular cloud
	vfreq=(freq(9500)-freq)*c/freq(9500)
	a5=(disk2_res)/1.665
	line_prof5=-exp(-(vfreq/a5)^2)/(SQRT(!pi)*a5)
	line_prof5=SHIFT(line_prof5,-9500)
	conv_spec=REAL_PART(FFT(FFT(spectrum_disk_d13)*FFT(line_prof5),1))
        spectrum_disk_d13=conv_spec*total(spectrum_disk_d13)/total(conv_spec) ; <--SPECTRUM of DISK!!!
	conv_spec=REAL_PART(FFT(FFT(spectrum_disk_d)*FFT(line_prof5),1))
	spectrum_disk_d=conv_spec*total(spectrum_disk_d)/total(conv_spec) ; <--SPECTRUM of DISK!!!
        conv_spec=REAL_PART(FFT(FFT(spectrum_disk_v21_d)*FFT(line_prof5),1))
        spectrum_disk_v21_d=conv_spec*total(spectrum_disk_v21_d)/total(conv_spec) ; <--SPECTRUM of DISK!!!





	wind_res=20.0
	vfreq=(freq(9500)-freq)*c/freq(9500)
	a4=(wind_res)/1.665
	line_prof2=-exp(-(vfreq/a4)^2)/(SQRT(!pi)*a4)
	line_prof2=SHIFT(line_prof2,-9500)
	conv_spec=REAL_PART(FFT(FFT(spectrum_wind)*FFT(line_prof2),1))
	spectrum_wind=conv_spec*total(spectrum_wind)/total(conv_spec) ; <--SPECTRUM of WIND!!!


	env_res=4.0
	vfreq=(freq(9500)-freq)*c/freq(9500)
	a7=(env_res)/1.665
	line_prof7=-exp(-(vfreq/a7)^2)/(SQRT(!pi)*a7)
	line_prof7=SHIFT(line_prof7,-9500)
	conv_spec=REAL_PART(FFT(FFT(spectrum_env)*FFT(line_prof7),1))
	spectrum_env=conv_spec*total(spectrum_env)/total(conv_spec) ; <--SPECTRUM of ENVELOPE!!!

	env2_res=5.0
	vfreq=(freq(9500)-freq)*c/freq(9500)
	a9=(env2_res)/1.665
	line_prof9=-exp(-(vfreq/a9)^2)/(SQRT(!pi)*a9)
	line_prof9=SHIFT(line_prof9,-9500)
	conv_spec=REAL_PART(FFT(FFT(spectrum_disk_cc)*FFT(line_prof9),1))
	spectrum_disk_cc=conv_spec*total(spectrum_disk_cc)/total(conv_spec) ; <--SPECTRUM of ENVELOPE!!!
        conv_spec=REAL_PART(FFT(FFT(spectrum_disk_cc13)*FFT(line_prof9),1))
        spectrum_disk_cc13=conv_spec*total(spectrum_disk_cc13)/total(conv_spec) ; <--SPECTRUM of ENVELOPE!!!
        conv_spec=REAL_PART(FFT(FFT(spectrum_disk_v21_cc)*FFT(line_prof9),1))
        spectrum_disk_v21_cc=conv_spec*total(spectrum_disk_v21_cc)/total(conv_spec) ; <--SPECTRUM of ENVELOPE!!!




	spectrum=spectrum_disk_v21_cc*spectrum_disk_v21_c*spectrum_disk_v21_d*spectrum_disk_cc*spectrum_disk_cc13*spectrum_env*spectrum_disk_c*spectrum_disk_c13*spectrum_disk_d*spectrum_disk_d13*spectrum_wind; + (1-exp(-tau_diskem12))*B_T12 + (1-exp(-tau_diskem13))*B_T13 + (1-exp(-tau_diskem18))*B_T18)

	;convolve the precise spectrum we calculate with the instrument resolution
	vfreq=(freq(9500)-freq)*c/freq(9500)
	a2=(inst_res)/1.665
	line_prof=-exp(-(vfreq/a2)^2)/(SQRT(!pi)*a2)
	line_prof=SHIFT(line_prof,-9500)
	conv_spec=REAL_PART(FFT(FFT(spectrum)*FFT(line_prof),1))
	conv_spec=conv_spec*total(spectrum)/total(conv_spec) ; <--FINAL SPECTRUM!!!


	;12CO chisq fits

	diff=ratio-INTERPOL(conv_spec,freq-44.15*freq/2.9979e5,f)
	count_tot=0
	;goto,step_4win
	;index=WHERE(f GE 2000 and f LE 2200 and finite(diff) EQ 1,count)
	index1=WHERE(f GE 1958.00 and f LE 1958.60 and finite(diff) EQ 1,count)
	chisq=TOTAL(diff(index1)^2/.04^2)
	count_tot=count_tot+count

	index2=WHERE(f GE 1978.40 and f LE 1978.85 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index2)^2/.04^2)
	count_tot=count_tot+count

	index3=WHERE(f GE 2008.0 and f LE 2008.50 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index3)^2/.04^2)
	count_tot=count_tot+count

	index4=WHERE(f GE 2031.80 and f LE 2032.35 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index4)^2/.04^2)
	count_tot=count_tot+count

	index5=WHERE(f GE 2068.2 and f LE 2068.85 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index5)^2/.04^2)
	count_tot=count_tot+count

	index6=WHERE(f GE 2085.7 and f LE 2086.30 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index6)^2/.04^2)
	count_tot=count_tot+count


	index7=WHERE(f GE 2106.80 and f LE 2107.35 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index7)^2/.04^2)
	count_tot=count_tot+count


	index8=WHERE(f GE 2135.0 and f LE 2135.5  and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index8)^2/.04^2)
	count_tot=count_tot+count


	rbig_model=conv_spec
	fbig_model=freq-44.15*freq/2.9979e5;fitpar12_10_avg(10)*freq/2.9979e5


	goto,step_13co

	cgps_open,'model_12CO_I.ps'
	cgdisplay,xsize=1100,ysize=750
	;goto,step_endR
	;wset,0
	str='P41'
	index=WHERE(id12_10(incl_line_12_10) EQ str)
	cgplot,fbig,ratio,xra=[f12_10(incl_line_12_10(index))+0.1,f12_10(incl_line_12_10(index))-0.7],xs=1,xtickinterval=0.4,yra=[-0.5,1.4],ys=1,ytitle='Normalized Flux',xtitle='Wavenumber (cm$\up-1$)',position=[0.10,0.59,0.31,0.99],charsize=1,psym=10
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index1),diff(index1)-0.2,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1

	str='P37'
	index=WHERE(id12_10(incl_line_12_10) EQ str)
	cgplot,fbig,rbig,xra=[f12_10(incl_line_12_10(index))+0.1,f12_10(incl_line_12_10(index))-0.7],xs=1,xtickinterval=0.4,ytickformat="(A1)",yra=[-0.5,1.4],ys=1,position=[0.31,0.59,0.52,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index2),diff(index2)-0.2,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1

	str='P31'
	index=WHERE(id12_10(incl_line_12_10) EQ str)
	cgplot,fbig,rbig,xra=[f12_10(incl_line_12_10(index))+0.1,f12_10(incl_line_12_10(index))-0.7],xs=1,xtickinterval=0.4,ytickformat="(A1)",yra=[-0.5,1.4],ys=1,position=[0.52,0.59,0.73,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index3),diff(index3)-0.2,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1

	str='P26'
	index=WHERE(id12_10(incl_line_12_10) EQ str)
	cgplot,fbig,rbig,xra=[f12_10(incl_line_12_10(index))+0.1,f12_10(incl_line_12_10(index))-0.7],xs=1,xtickinterval=0.4,ytickformat="(A1)",yra=[-0.5,1.4],ys=1,position=[0.73,0.59,0.94,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index4),diff(index4)-0.2,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1

	str='P18'
	index=WHERE(id12_10(incl_line_12_10) EQ str)
	cgplot,fbig,rbig,xra=[f12_10(incl_line_12_10(index))+0.1,f12_10(incl_line_12_10(index))-0.7],xs=1,xtickinterval=0.4,yra=[-0.5,1.4],ys=1,ytitle='Normalized Flux',position=[0.10,0.10,0.31,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index5),diff(index5)-0.2,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1

	str='P14'
	index=WHERE(id12_10(incl_line_12_10) EQ str)
	cgplot,fbig,rbig,xra=[f12_10(incl_line_12_10(index))+0.1,f12_10(incl_line_12_10(index))-0.7],xs=1,xtickinterval=0.3,ytickformat="(A1)",yra=[-0.5,1.4],ys=1,position=[0.31,0.10,0.52,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index6),diff(index6)-0.2,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1




	str='P9'
	index=WHERE(id12_10(incl_line_12_10) EQ str)
	cgplot,fbig,rbig,xra=[f12_10(incl_line_12_10(index))+0.1,f12_10(incl_line_12_10(index))-0.7],xs=1,xtickinterval=0.4,ytickformat="(A1)",yra=[-0.5,1.4],ys=1,position=[0.52,0.10,0.73,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index7),diff(index7)-0.2,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1


	str='P2'
	index=WHERE(id12_10(incl_line_12_10) EQ str)
	cgplot,fbig,rbig,xra=[f12_10(incl_line_12_10(index))+0.1,f12_10(incl_line_12_10(index))-0.7],xs=1,xtickinterval=0.4,ytickformat="(A1)",yra=[-0.5,1.4],ys=1,position=[0.73,0.10,0.94,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index8),diff(index8)-0.2,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1
	;goto,step_end

	cgps_close,/PNG



	;cgps_open,'missingcomp.ps'
	;cgdisplay,xsize=600,ysize=1000
	;str='P31'
        ;index=WHERE(id12_10(incl_line_12_10) EQ str)
	;vbig1=fbig
	;FOR i=0,N_ELEMENTS(fbig)-1 DO vbig1(i)=(f12_10(incl_line_12_10(index)) - fbig(i))*2.9979e5/f12_10(incl_line_12_10(index))
        ;cgplot,vbig1-44.15+17.3,rbig,xra=[-40,70],xs=1,xtickinterval=10.0,yra=[-0.5,1.4],ys=1,ytitle='Normalized Flux',position=[0.20,0.59,0.94,0.99],charsize=1,psym=10,/NOERASE
	;vbig1_model=fbig_model
	;FOR i=0,N_ELEMENTS(fbig_model)-1 DO vbig1_model(i)=(f12_10(incl_line_12_10(index)) - fbig_model(i))*2.9979e5/f12_10(incl_line_12_10(index))
	;cgoplot,vbig1_model-44.15+17.3,rbig_model,color='goldenrod',/OVERPLOT
	;cgoplot,[17.3,17.3],[-2,2],linestyle=2,/OVERPLOT
        ;cgoplot,vbig1(index3)-44.15+17.3,diff(index3),color='goldenrod',psym=10,/OVERPLOT
	;cgtext,-30,1.2,str,charsize=1

	;str='P2'
        ;index=WHERE(id12_10(incl_line_12_10) EQ str)
        ;vbig1=fbig
        ;FOR i=0,N_ELEMENTS(fbig)-1 DO vbig1(i)=(f12_10(incl_line_12_10(index)) - fbig(i))*2.9979e5/f12_10(incl_line_12_10(index))
        ;cgplot,vbig1-44.15+17.3,rbig,xra=[-40,70],xs=1,xtickinterval=10.0,yra=[-0.5,1.4],ys=1,ytitle='Normalized Flux',xtitle='Heliocentric Velocity (km s$\up-1$)',position=[0.20,0.13,0.94,0.53],charsize=1,psym=10,/NOERASE
        ;vbig1_model=fbig_model
        ;FOR i=0,N_ELEMENTS(fbig_model)-1 DO vbig1_model(i)=(f12_10(incl_line_12_10(index)) - fbig_model(i))*2.9979e5/f12_10(incl_line_12_10(index))
        ;cgoplot,vbig1_model-44.15+17.3,rbig_model,color='goldenrod',/OVERPLOT
        ;cgoplot,[17.3,17.3],[-2,2],linestyle=2,/OVERPLOT
        ;cgoplot,vbig1(index8)-44.15+17.3,diff(index8),color='goldenrod',psym=10,/OVERPLOT
        ;cgtext,-30,1.2,str,charsize=1

	
	;cgplot,fbig,rbig,xra=[f12_10(incl_line_12_10(index))-0.7,f12_10(incl_line_12_10(index))+0.1],xs=1,xtickinterval=0.3,yra=[-0.5,1.4],ys=1,xtitle='Wavenumber (cm$\up-1$)',ytitle='Normalized Flux',position=[0.20,0.13,0.94,0.53],charsize=1,psym=10,/NOERASE
        ;cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
        ;cgoplot,fbig(index8),diff(index8),color='goldenrod',psym=10,/OVERPLOT
        ;cgtext,f12_10(incl_line_12_10(index))-0.1,1.2,str,charsize=1
	;pause
	;cgps_close,/PNG


	step_13co:
	goto,step_c18o
	;13CO chisq fits
	index1=WHERE(f GE 1994.10 and f LE 1994.50 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index1)^2/.04^2)
	count_tot=count_tot+count

	index2=WHERE(f GE 2011.60 and f LE 2012.00 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index2)^2/.04^2)
	count_tot=count_tot+count

	index3=WHERE(f GE 2032.8 and f LE 2033.30 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index3)^2/.04^2)
	count_tot=count_tot+count

	index4=WHERE(f GE 2057.30 and f LE 2057.60 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index4)^2/.04^2)
	count_tot=count_tot+count

	index5=WHERE(f GE 2120.3 and f LE 2120.8 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index5)^2/.04^2)
	count_tot=count_tot+count

	index6=WHERE(f GE 2113.4 and f LE 2113.8  and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index6)^2/.04^2)
	count_tot=count_tot+count


	index7=WHERE(f GE 2084.40 and f LE 2084.80 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index7)^2/.04^2)
	count_tot=count_tot+count


	index8=WHERE(f GE 2106.4 and f LE 2106.7  and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index8)^2/.04^2)
	count_tot=count_tot+count
	
	;goto,stepII_chisq

	;cgps_open,'model_13CO.ps'
	cgdisplay,xsize=1100,ysize=750
	;wset,1
	str='P25'
	index=WHERE(id13_10(incl_line_13_10) EQ str)
	cgplot,fbig,ratio,xra=[f13_10(incl_line_13_10(index))+0.1,f13_10(incl_line_13_10(index))-0.7],xs=1,yra=[-0.5,1.4],ys=1,xtickinterval=0.4,ytitle='Normalized Flux',xtitle='Wavenumber (cm$\up-1$)',position=[0.10,0.59,0.31,0.99],charsize=1,psym=10
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index1),diff(index1)-0.1,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f13_10(incl_line_13_10(index))-0.1,1.2,str,charsize=1

	str='P21'
	index=WHERE(id13_10(incl_line_13_10) EQ str)
	cgplot,fbig,rbig,xra=[f13_10(incl_line_13_10(index))+0.1,f13_10(incl_line_13_10(index))-0.7],xs=1,yra=[-0.5,1.4],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.31,0.59,0.52,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index2),diff(index2)-0.1,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f13_10(incl_line_13_10(index))-0.1,1.2,str,charsize=1

	str='P16'
	index=WHERE(id13_10(incl_line_13_10) EQ str)
	cgplot,fbig,rbig,xra=[f13_10(incl_line_13_10(index))+0.1,f13_10(incl_line_13_10(index))-0.7],xs=1,yra=[-0.5,1.4],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.52,0.59,0.73,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index3),diff(index3)-0.1,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f13_10(incl_line_13_10(index))-0.1,1.2,str,charsize=1

	str='P10'
	index=WHERE(id13_10(incl_line_13_10) EQ str)
	cgplot,fbig,rbig,xra=[f13_10(incl_line_13_10(index))+0.1,f13_10(incl_line_13_10(index))-0.7],xs=1,yra=[-0.5,1.4],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.73,0.59,0.94,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index4),diff(index4)-0.1,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f13_10(incl_line_13_10(index))-0.1,1.2,str,charsize=1

	chisq_test=chisq/count_tot

	str='R6'
	index=WHERE(id13_10(incl_line_13_10) EQ str)
	cgplot,fbig,rbig,xra=[f13_10(incl_line_13_10(index))+0.1,f13_10(incl_line_13_10(index))-0.7],xs=1,yra=[-0.5,1.4],ys=1,xtickinterval=0.4,ytitle='Normalized Flux',position=[0.10,0.10,0.31,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index5),diff(index5)-0.1,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f13_10(incl_line_13_10(index))-0.1,1.2,str,charsize=1


	str='R4'
	index=WHERE(id13_10(incl_line_13_10) EQ str)
	cgplot,fbig,rbig,xra=[f13_10(incl_line_13_10(index))+0.1,f13_10(incl_line_13_10(index))-0.7],xs=1,yra=[-0.5,1.4],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.31,0.10,0.52,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index6),diff(index6)-0.1,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f13_10(incl_line_13_10(index))-0.1,1.2,str,charsize=1




	str='P3'
	index=WHERE(id13_10(incl_line_13_10) EQ str)
	cgplot,fbig,rbig,xra=[f13_10(incl_line_13_10(index))+0.1,f13_10(incl_line_13_10(index))-0.7],xs=1,yra=[-0.5,1.4],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.52,0.10,0.73,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index7),diff(index7)-0.1,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f13_10(incl_line_13_10(index))-0.1,1.2,str,charsize=1


	str='R2'
	index=WHERE(id13_10(incl_line_13_10) EQ str)
	cgplot,fbig,rbig,xra=[f13_10(incl_line_13_10(index))+0.1,f13_10(incl_line_13_10(index))-0.7],xs=1,yra=[-0.5,1.4],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.73,0.10,0.94,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index8),diff(index8)-0.1,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f13_10(incl_line_13_10(index))-0.1,1.2,str,charsize=1
	;cgps_close,/PNG
	
	
	step_c18o:
	goto,step_12cov21
	;C18O chisq fits
	index1=WHERE(f GE 2145.65 and f LE 2146.00 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index1)^2/.04^2)
	count_tot=count_tot+count

	index2=WHERE(f GE 2132.85 and f LE 2133.30 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index2)^2/.04^2)
	count_tot=count_tot+count

	index3=WHERE(f GE 2057.6 and f LE 2057.85 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index3)^2/.04^2)
	count_tot=count_tot+count

	index4=WHERE(f GE 2126.40 and f LE 2126.70 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index4)^2/.04^2)
	count_tot=count_tot+count

	index5=WHERE(f GE 2119.70 and f LE 2120.10 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index5)^2/.04^2)
	count_tot=count_tot+count

	index6=WHERE(f GE 2109.50 and f LE 2109.75  and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index6)^2/.04^2)
	count_tot=count_tot+count


	index7=WHERE(f GE 2105.95 and f LE 2106.30 and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index7)^2/.04^2)
	count_tot=count_tot+count


	index8=WHERE(f GE 2084.3 and f LE 2084.5  and finite(diff) EQ 1,count)
	chisq=chisq + TOTAL(diff(index8)^2/.04^2)
	count_tot=count_tot+count
	chisq=chisq/count_tot
	print,'CHISQ: ' + STRTRIM(STRING(chisq),1)
	
	




	;goto,step_end

	cgps_open,'model_C18O.ps'
	cgdisplay,xsize=1100,ysize=750
	;wset,2
	str='R15'
	index=WHERE(id18_10(incl_line_18_10) EQ str)
	cgplot,fbig,ratio,xra=[f18_10(incl_line_18_10(index))+0.1,f18_10(incl_line_18_10(index))-0.7],xs=1,yra=[0.4,1.3],ys=1,xtitle='Wavenumber (cm$\up-1$)',ytitle='Normalized Flux',xtickinterval=0.4,position=[0.10,0.59,0.31,0.99],charsize=1,psym=10
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index1),diff(index1)+0.6,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f18_10(incl_line_18_10(index))-0.1,1.2,str,charsize=1

	str='R11'
	index=WHERE(id18_10(incl_line_18_10) EQ str)
	cgplot,fbig,rbig,xra=[f18_10(incl_line_18_10(index))+0.1,f18_10(incl_line_18_10(index))-0.7],xs=1,yra=[0.4,1.3],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.31,0.59,0.52,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index2),diff(index2)+0.6,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f18_10(incl_line_18_10(index))-0.1,1.2,str,charsize=1

	str='P9'
	index=WHERE(id18_10(incl_line_18_10) EQ str)
	cgplot,fbig,rbig,xra=[f18_10(incl_line_18_10(index))+0.1,f18_10(incl_line_18_10(index))-0.7],xs=1,yra=[0.4,1.3],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.52,0.59,0.73,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index3),diff(index3)+0.6,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f18_10(incl_line_18_10(index))-0.1,1.2,str,charsize=1

	str='R9'
	index=WHERE(id18_10(incl_line_18_10) EQ str)
	cgplot,fbig,rbig,xra=[f18_10(incl_line_18_10(index))+0.1,f18_10(incl_line_18_10(index))-0.7],xs=1,yra=[0.4,1.3],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.73,0.59,0.94,0.99],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index4),diff(index4)+0.6,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f18_10(incl_line_18_10(index))-0.1,1.2,str,charsize=1

	str='R7'
	index=WHERE(id18_10(incl_line_18_10) EQ str)
	cgplot,fbig,rbig,xra=[f18_10(incl_line_18_10(index))+0.1,f18_10(incl_line_18_10(index))-0.7],xs=1,yra=[0.4,1.3],ys=1,xtickinterval=0.4,ytitle='Normalized Flux',position=[0.10,0.10,0.31,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index5),diff(index5)+0.6,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f18_10(incl_line_18_10(index))-0.1,1.2,str,charsize=1

	str='R4'
	index=WHERE(id18_10(incl_line_18_10) EQ str)
	cgplot,fbig,rbig,xra=[f18_10(incl_line_18_10(index))+0.1,f18_10(incl_line_18_10(index))-0.7],xs=1,yra=[0.4,1.3],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.31,0.10,0.52,0.50],charsize=1,psym=10,/NOERASE
	cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
	cgoplot,fbig(index6),diff(index6)+0.6,color='goldenrod',psym=10,/OVERPLOT
	cgtext,f18_10(incl_line_18_10(index))-0.1,1.2,str,charsize=1

	str='R3'
        index=WHERE(id18_10(incl_line_18_10) EQ str)
        cgplot,fbig,rbig,xra=[f18_10(incl_line_18_10(index))+0.1,f18_10(incl_line_18_10(index))-0.7],xs=1,yra=[0.4,1.3],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.52,0.10,0.73,0.50],charsize=1,psym=10,/NOERASE
        cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
        cgoplot,fbig(index7),diff(index7)+0.6,color='goldenrod',psym=10,/OVERPLOT
        cgtext,f18_10(incl_line_18_10(index))-0.1,1.2,str,charsize=1

        str='P2'
        index=WHERE(id18_10(incl_line_18_10) EQ str)
        cgplot,fbig,rbig,xra=[f18_10(incl_line_18_10(index))+0.1,f18_10(incl_line_18_10(index))-0.7],xs=1,yra=[0.4,1.3],ys=1,xtickinterval=0.4,ytickformat="(A1)",position=[0.73,0.10,0.94,0.50],charsize=1,psym=10,/NOERASE
        cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
        cgoplot,fbig(index8),diff(index8)+0.6,color='goldenrod',psym=10,/OVERPLOT
        cgtext,f18_10(incl_line_18_10(index))-0.1,1.2,str,charsize=1



	cgps_close,/PNG
;ENDFOR
;goto,step_endII

step_12cov21:
;12CO v2-1

index1=WHERE(f GE 2178.60 and f LE 2179.00 and finite(diff) EQ 1,count)
chisq=chisq + TOTAL(diff(index1)^2/.04^2)
count_tot=count_tot+count

index2=WHERE(f GE 2162.50 and f LE 2162.85 and finite(diff) EQ 1,count)
chisq=chisq + TOTAL(diff(index2)^2/.04^2)
count_tot=count_tot+count

index3=WHERE(f GE 2138.35 and f LE 2138.70 and finite(diff) EQ 1,count)
chisq=chisq + TOTAL(diff(index3)^2/.04^2)
count_tot=count_tot+count

index4=WHERE(f GE 2096.80 and f LE 2097.20 and finite(diff) EQ 1,count)
chisq=chisq + TOTAL(diff(index4)^2/.04^2)
count_tot=count_tot+count

index5=WHERE(f GE 2088.8 and f LE 2089.20 and finite(diff) EQ 1,count)
chisq=chisq + TOTAL(diff(index5)^2/.04^2)
count_tot=count_tot+count

index6=WHERE(f GE 2084.8 and f LE 2085.2  and finite(diff) EQ 1,count)
chisq=chisq + TOTAL(diff(index6)^2/.04^2)
count_tot=count_tot+count


index7=WHERE(f GE 2055.50 and f LE 2055.95 and finite(diff) EQ 1,count)
chisq=chisq + TOTAL(diff(index7)^2/.04^2)
count_tot=count_tot+count


index8=WHERE(f GE 2051.2 and f LE 2051.60 and finite(diff) EQ 1,count)
chisq=chisq + TOTAL(diff(index8)^2/.04^2)
count_tot=count_tot+count

cgps_open,'model_12COv21.ps'
cgdisplay,xsize=1100,ysize=750


RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/CO_molecdat.dat'
f12_21 = X12CO21(6,*)




;wset,3
str='R17'
index=WHERE(id12_21(incl_line_12_21) EQ str)
;print,index
print,f12_21(incl_line_12_21(index))+0.1
cgplot,fbig,ratio,xra=[f12_21(incl_line_12_21(index))+0.1,f12_21(incl_line_12_21(index))-0.7],xs=1,yra=[-0.2,1.3],ys=1,xtickinterval=0.3,xtitle='Wavenumber (cm$\up-1$)',ytitle='Normalized Flux',position=[0.10,0.59,0.31,0.99],charsize=1,psym=10
cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
cgoplot,fbig(index1),diff(index1),color='goldenrod',/OVERPLOT
cgtext,f12_21(incl_line_12_21(index))-0.1,1.2,str,charsize=1

str='R12'
index=WHERE(id12_21(incl_line_12_21) EQ str)
;print,index
print,f12_21(incl_line_12_21(index))+0.1
cgplot,fbig,rbig,xra=[f12_21(incl_line_12_21(index))+0.1,f12_21(incl_line_12_21(index))-0.7],xs=1,yra=[-0.2,1.3],ys=1,xtickinterval=0.3,ytickformat="(A1)",position=[0.31,0.59,0.52,0.99],charsize=1,psym=10,/NOERASE
cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
cgoplot,fbig(index2),diff(index2),color='goldenrod',/OVERPLOT
cgtext,f12_21(incl_line_12_21(index))-0.1,1.2,str,charsize=1

str='R5'
index=WHERE(id12_21(incl_line_12_21) EQ str)
cgplot,fbig,rbig,xra=[f12_21(incl_line_12_21(index))+0.1,f12_21(incl_line_12_21(index))-0.7],xs=1,yra=[-0.2,1.3],ys=1,xtickinterval=0.3,ytickformat="(A1)",position=[0.52,0.59,0.73,0.99],charsize=1,psym=10,/NOERASE
cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
cgoplot,fbig(index3),diff(index3),color='goldenrod',/OVERPLOT
cgtext,f12_21(incl_line_12_21(index))-0.1,1.2,str,charsize=1




str='P5'
index=WHERE(id12_21(incl_line_12_21) EQ str)
cgplot,fbig,rbig,xra=[f12_21(incl_line_12_21(index))+0.1,f12_21(incl_line_12_21(index))-0.7],xs=1,yra=[-0.2,1.3],ys=1,xtickinterval=0.3,ytickformat="(A1)",position=[0.73,0.59,0.94,0.99],charsize=1,psym=10,/NOERASE
cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
cgoplot,fbig(index4),diff(index4),color='goldenrod',/OVERPLOT
cgtext,f12_21(incl_line_12_21(index))-0.1,1.2,str,charsize=1

str='P7'
index=WHERE(id12_21(incl_line_12_21) EQ str)
cgplot,fbig,rbig,xra=[f12_21(incl_line_12_21(index))+0.1,f12_21(incl_line_12_21(index))-0.7],xs=1,yra=[-0.2,1.3],ys=1,ytitle='Normalized Flux',xtickinterval=0.3,position=[0.10,0.10,0.31,0.50],charsize=1,psym=10,/NOERASE
cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
cgoplot,fbig(index5),diff(index5),color='goldenrod',/OVERPLOT
cgtext,f12_21(incl_line_12_21(index))-0.1,1.2,str,charsize=1

str='P8'
index=WHERE(id12_21(incl_line_12_21) EQ str)
cgplot,fbig,rbig,xra=[f12_21(incl_line_12_21(index))+0.1,f12_21(incl_line_12_21(index))-0.7],xs=1,yra=[-0.2,1.3],ys=1,xtickinterval=0.3,ytickformat="(A1)",position=[0.31,0.10,0.52,0.50],charsize=1,psym=10,/NOERASE
cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
cgoplot,fbig(index6),diff(index6),color='goldenrod',/OVERPLOT
cgtext,f12_21(incl_line_12_21(index))-0.1,1.2,str,charsize=1




str='P15'
index=WHERE(id12_21(incl_line_12_21) EQ str)
cgplot,fbig,rbig,xra=[f12_21(incl_line_12_21(index))+0.1,f12_21(incl_line_12_21(index))-0.7],xs=1,yra=[-0.2,1.3],ys=1,xtickinterval=0.3,ytickformat="(A1)",position=[0.52,0.10,0.73,0.50],charsize=1,psym=10,/NOERASE
cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
cgoplot,fbig(index7),diff(index7),color='goldenrod',/OVERPLOT
cgtext,f12_21(incl_line_12_21(index))-0.1,1.2,str,charsize=1


str='P16'
index=WHERE(id12_21(incl_line_12_21) EQ str)
cgplot,fbig,rbig,xra=[f12_21(incl_line_12_21(index))+0.1,f12_21(incl_line_12_21(index))-0.7],xs=1,yra=[-0.2,1.3],ys=1,xtickinterval=0.3,ytickformat="(A1)",position=[0.73,0.10,0.94,0.50],charsize=1,psym=10,/NOERASE
cgoplot,fbig_model,rbig_model,color='goldenrod',/OVERPLOT
cgoplot,fbig(index8),diff(index8),color='goldenrod',/OVERPLOT
cgtext,f12_21(incl_line_12_21(index))-0.1,1.2,str,charsize=1


cgps_close,/PNG
;pause

step_end:

;tmp=[chisq,veil_a,veil_cc,veil_cc13,veil_c,veil_d,ds12_a,ds12_b,ds12_c,ds12_cc,ds12_d,ds13_b,ds13_c,ds13_cc,ds13_d,ds18_c,ds18_cc,ds18_d,trot12a,trot12b,trot12c,trot12cc,trot12d,N12a,N12b,N12c,N12cc,N12d,int_wth12a,int_wth12b,int_wth12c,int_wth12cc,int_wth12d,trot1221,N1221,int_wth1221,trot13b,trot13c,trot13cc,trot13d,N13b,N13c,N13cc,N13d,int_wth13b,int_wth13c,int_wth13cc,int_wth13d,trot18c,trot18cc,trot18d,N18c,N18cc,N18d,int_wth18c,int_wth18cc,int_wth18d]

;print, chisq

;find decent fits
;IF chisq LT 5.0 THEN BEGIN
;OPENW,1,'par_gvtau_trot12cc.txt',/append
;PRINTF,1,tmp,FORMAT='(F11.6,5F6.2,12F11.4,5F11.4,5E13.4,5F11.6,F11.4,E13.4,F11.6,4F11.4,4E13.4,4F11.6,3F11.4,3E13.4,3F11.6)'
;CLOSE,1
;ENDIF


ENDFOR

goto,step_endII







step_endII:


END
