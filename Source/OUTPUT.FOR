c pot_ET is used to save the accumulated potential transpiration in the previous N steps
c p_ET is used to save the difference between the accumulated transpiration in the 
c previous N step and N+1 step; thus in fact it is the potential transpiration rate per time step
      subroutine Output ()
      Include 'public.ins'
      Include 'PuWeath.ins'
	  Include 'PuPlant.ins'

      Integer Plevel,e
      
      Dimension temp(1000)
      Real PotentialSoilEvaporation, WidthT,Period
      Real SeasonRainfall, SeasonPotenTranspiration,
     &     SeasonActualTranspiration, SeasonPotenSoilEvap,
     &     SeasonActualSoilEvap,SeasonActualInfil
      Character*10 date, Stime
      Common /Outcom/NumPoint(2000),NodNum(2000,60),
     &        NumCell(1000),NumElPrint(1000,60),ModNum,
     &        Period, ETO, ActualTrans1,ActualSoilEvap,
     &        pot_ET,pot_ET_EO,bottomFlux, ActualInfil,timer,
     &        FluxAct,Runoff,PotentialSoilEvaporation,
     &        CumRain,Runoff_Terminal,
     &        SeasonRainfall, SeasonPotenTranspiration,
     &        SeasonActualTranspiration, SeasonPotenSoilEvap,
     &        SeasonActualSoilEvap,SeasonActualInfil,CumulativeFlux

      t=Time
      
      If(lInput.eq.1) then
C
C  Initialize 
C
        If (DailyOutput.eq.1) Period=1.0
        If (HourlyOutput.eq.1) Period=1.0/(24.0)
        NumMod=NumMod+1
        ModNum=NumMod
        tNext(ModNum) = time

C*** FROM HourlyOut.for BEGIN      
        timer=0
        bottomFlux=0.0
        fluxAct=0.0
        Runoff=0.0
        Runoff_Terminal=0.0
        ActualSoilEvap=0
        ActualInfil=0
        PotentialSoilEvaporation=0.0
        CumRain=0.0
        SeasonRainfall =0.0
        SeasonActualInfil=0.0
        SeasonPotenTranspiration=0.0
        SeasonActualTranspiration=0.0
        SeasonPotenSoilEvap=0.0
        SeasonActualSoilEvap=0.0
	  ActualTrans1=0.0
	  pot_ET=0.0
	  pot_ET_EO=0.0
        
	  Open(81,file=NodeGraphics)
	  open(82,file=ElemGraphics)  
	  Open(83,file=SurfaceGraphics)
        Open(84,file=FluxGraphics)
         
      Write(81,5) "Date_time", "Date", "X", "Y", "hNew","thNew", "VZ", 
     &   "VX","Q", "ConcN", "Temp","GasConc" ,"Capa", "Cond" 
     
      Write(82,6) "Date_time", "Date", "X", "Y", "Elem", "RTWT",
     &         "TotalRT", "WaterSink", "NitSink", "GasSink", "Area"
         
        Write(83,20) "Date_time", "Date", "PSoilEvap", 
     &         "ASoilEVap", "PET_PEN",
     &         "PE_T_int","transp","CumRain","infil","FLuxAct",
     &         "BotFlux", "Runoff","cover","PSIM", "SeasPSoEv",
     &         "SeasASoEv",
     &        "SeasPTran","SeasATran","SeasRain","SeasInfil" 
     
       Write(84,20) "time","Date","Flux", "ActualSoilEvap",
     &                 "CumulativeFlux"
       
       CumulativeFlux=0
       Endif     !end initialization
C
C Routine calculations every time step
C

C Potential values:

*  add only if surface is not shaded
* note that VarBW(i,3) is the adjusted et or rainfall - i.e., sum of rain + et, not separated
* we may have to consider changing this in the future
* 6/2016 DT changed CumRain sum to use varBW(i,1) which is actual rain.
* Actual infiltration, after accounting for evap is calculated in ActualInfil
          Do i=1,NumBp
	          n=KXB(i)
             If (VarBW(i,3).GE.0) Then
      	        PotentialSoilEvaporation=PotentialSoilEvaporation+
     & VarBW(i,3)*width(i)*step
	         EndIf

	       End Do
             
c CCCCCCCCCCCCCCCCCCCC We change the VarBW as Rain+Runoff
c CCCCCCCCCCCCCCCCCCCC We have to record the Rain after the each run of hourwea
c Since we edit the VarBW as Rain+Runoff(i-1), this should be cumu infiltration or so
c                CumRain=CumRain+varBW(i,1)*width(i)*step
c This is the real cumuRain
          do i=1,NumBp
             n=KXB(i)
             if (codeW(n).eq.-4) Then
              CumRain=CumRain+VarBW_Old_1(i)*width(i)*step
             endif
          enddo
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	
C Actual values

CDT note that flow takes place in seepage faces when the boundary code is positive - this indicates a fixed
C  pressure head and bottomFlux is calculated CodeW >0 indicates head is set at 0 and flux is calculated.
         timer=timer+step
      	 do i=1,NumBp
	     n=KXB(i)
	     
	     
	     if ((CodeW(n).eq.(-7)).or.(CodeW(n).eq.(2))
     &             .or.(CodeW(n).eq.(1))) then
	        bottomFlux=bottomFlux+Q(n)*step
	        if (bottomFlux.gt.0) then 
	          iii=1
	          endif
	        else
	         if (CodeW(n).EQ.-2) then
	       	         iii=1
	       	  endif
           endif
c case for evaporation
           if ((abs(CodeW(n)).eq.4).and.(VarBw(i,3).GT.0)) then
              ActualSoilEvap=ActualSoilEvap-min(0.0,Q(n))*step
            endif
c case for rainfall
           if (abs(CodeW(n)).eq.4)
     !            then
              ActualInfil=ActualInfil+max(0.0,Q(n))*step
              Runoff=Runoff+RO(n)*step
c              RO(n)=0.0   ! we will use this to store ponded water for next time period
              FluxAct=FluxAct+max(0.0,QAct(n))*step
        
            endif

	  EndDo

c CCCCCCCCCCCCC Zhuangji: Need to collect the runoff from point 7
c CCCCCCCCCCCCC point 7 is manually determined, but need to change if the runoff point is given.
c              Runoff_Terminal=Runoff_Terminal+RO(7)*step
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  
	
C*** FROM HourlyOut.for END

11    If(abs(time-tNext(ModNum)).lt.0.001*Step.OR.lInput.eq.1) then
          iday=int(t)
        call caldat(iday,mm,id,iyyy) 
        write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
          Call ArrNodOut (Date,x,y,NumNP,NumEl,t,NumElD,hNew,
     !                   ThNew,vz,vx,Q,Conc(1,1),Tmpr,g(1,1),Cap,Con)

        If(NShoot.ne.0) then
             Call ArrElemOut (Date,x,y,NumNP,NumEl,t,KX,NumElD,RTWT,
     !                      RUTDEN,Sink,cSink(1,1),gSink(1,1),YRL,Area,
     !                      MatNumE)
            !Call OutToDB(Sink(1))
        Endif
        tNext(ModNum) = time + Period/10.0  ! every 10 minutes
C      
* This code operates hourly to output data on the surface nodes
	  
	  WidthT=0.0
        Do i= 1, NumBP
           if (abs(CodeW(KXB(i))).eq.4) then
           	  WidthT=WidthT+width(i)
           	 endif
        Enddo


        
c CCCCCCCCCCCCCCCCCCC Zhuangji's Surface calculation CCCCCCCCCCCCCCCC
          Surface_Flux=Surface_Flux+FluxAct
          Bottom_Flux=Bottom_Flux+bottomFlux
          Runoff_Flux=Runoff_Flux+Runoff
          Rain_Flux=Rain_Flux+CumRain
c we consider the runoff from the terminal edge of the soil profile
c          Runoff_Terminal_Flux=Runoff_Terminal_Flux+Runoff_Terminal
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

* Anything calculated from VarBW should not be divided by 24 
* output as mm per plant, essentially dividing by popslab
* VarBW is cm/cm2 per day. Any fluxes from the finite element code 
*  are cm of water in the width of the nodes per day. Thus these fluxes are not
* multiplied by width

              PotentialSoilEvaporation=PotentialSoilEvaporation
     &        /EOMult/poprow*100
              ActualSoilEvap=ActualSoilEvap/EOMult/poprow*100
              !ActualInfil=ActualInfil/EOMult/poprow*100
              !FluxAct=FluxAct/EOMult/poprow*100
              !Runoff=Runoff/EOMult/poprow*100
              !CumRain=CumRain/EOMult/poprow*100
              !bottomFlux=bottomFlux/EOMult/poprow*100
              !(/RowSP/2) gets us per cm2 now for 10 minute period
              !*6 gives us one hour
              ActualInfil=ActualInfil/(RowSP/2)*10*24
              FluxAct=FluxAct/(RowSP/2)*10*24
              Runoff=Runoff/(RowSP/2)*10*24
              CumRain=CumRain/(RowSP/2)*10*24
              bottomFlux=bottomFlux/(RowSP/2)*10*24    !cm/cm2/day
              CumulativeFlux=CumulativeFlux+bottomFlux


CDT distribute runoff back into Q This makes it uniform. implement this later

!c         do n=1, NumBP
!c           ROForInfil(n)=Runoff*Width(n)
!c          end do
          
* dt - value of awupss and pot_et are per half slab as calclated
* in the water uptake code. EOMULT accounts for whether the slab
* is a full or half (plant in the middle - full or edge -half)
* multiply by 1/poprow * 100  to multiply be the amount of space between 
* plants (cm) the scales the slab to the area corresponding to one
* plant - hence the output is per plant
* 

       
           CurrentTrans=(AWUPSS-ActualTrans1)/EOMult/poprow*100
           p_ET=(potential_T-pot_ET)/EOMult/poprow*100     ! already per plant
           
           p_ET_EO=potential_T_EO-pot_ET_EO     
c potential_T is calculated as a summation of EOR*Step so it is per slab and is cumulative over the simulation
C EOR is calculated from ET demand in the crop model
c  pot_ET is the value of potential_T from the previous hour
           
C EO is the potential transpiration rate from Penman adjusted for plant cover, calculated in the weather subroutine
c  potential_T_EO is the summation of EO which is soil surface evap from EO and the Penman equation.
c pot_ET_EO is the EO summation from the previous time step.           
c EOR_p is a rough estimate of potential ET. p_ET is a better estimate using time step to 
c  integrate rather than use an average
           EOR_p=EOR*period/EOMult/poprow*100
 

             
C 
C    Accumulate seasonal amounts       
C 
          SeasonRainfall =SeasonRainfall+CumRain
          SeasonActualInfil=SeasonActualInfil+ActualInfil
          SeasonPotenTranspiration=SeasonPotenTranspiration+
     &       p_ET
          SeasonActualTranspiration=SeasonActualTranspiration+
     &       CurrentTrans
          SeasonPotenSoilEvap=SeasonPotenSoilEvap+
     &       PotentialSoilEvaporation
          SeasonActualSoilEvap=SeasonActualSoilEvap+ActualSoilEvap
          
          
         Write(83,10) t,Date,PotentialSoilEvaporation, 
     &               ActualSoilEvap,p_ET_EO,
     &              p_ET, CurrentTrans,CumRain,ActualInfil,FLuxAct,
     &              bottomFlux, Runoff,cover, PSIM, SeasonPotenSoilEvap,
     &              SeasonActualSoilEvap, SeasonPotenTranspiration,
     &               SeasonActualTranspiration,  SeasonRainfall,  
     &              SeasonActualInfil
     
        Write(84,10) t,date,bottomFlux,ActualSoilEvap,CumulativeFlux

	  bottomFlux=0
	  ActualSoilEvap=0
	  ActualInfil=0
	  ActualTrans1=AWUPSS
	  pot_ET=potential_T
	  pot_ET_EO=Potential_T_EO
	  timer=0
	  FluxAct=0.0
	  Runoff=0.0
        Runoff_Terminal=0.0
	  CumRain=0.0
	  PotentialSoilEvaporation=0.0
        
c CCCCCCCCCCCCCCCCCCCCCC Zhuangji's Water Balance Term CCCCCCCCCCCCCCCCCCCC
c        Surface_Flux=0.0
c        Bottom_Flux=0.0
c        Runoff_Flux=0.0
c        Soil_Water_Diff=0.0
c        Balance_Z=0.0
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        
	  Endif  !end loop for this timestep
C*** FROM HourlyOut.for END
           
      Return
C10    Call errmes(im,il)
  5   Format (1x,A16, T20,A4,T36,A1,T46,A1, T52,A4,T63,A5,T75,A2,T87,A2,
     &     T97,A1, T108, A5, T119, A4, T129, A7,T97,A1,T97,A1)   ! format for nodal data header
  6   Format (1x,A16,T20,A4,T36,A1,T46,A1,T56,A4,T62,A4,T70,A7, T81,A9,
     &        T93, A7, T105, A7, T115, A4)
 10   Format (1x,F17.6,A15,20F16.3)  !format for flux data
 20   Format (1x,20A16)              !format for flux header
      End        
C/***************************************************/
      Subroutine ArrElemOut(Date,x,y,NumNP,NumEl,t,KX,NumElD,
     !          RTWT,RUTDEN,SinkArr,cSinkArr,gSinkArr,YRL,Area, MatNumE)
      Character *10 Date
      Dimension RTWT(NumEl),Area(NumEl),x(NumNP),y(NumNP),
     !        KX(NumElD,4),RUTDEN(NumEl),YRL(NumEl),MatNumE(NumEl),
     !         SinkArr(NumEl),cSinkArr(NumEl),gSinkArr(NumEl)
        

        Do m=1,NumEl
          NUS=4
          If(KX(m,3).eq.KX(m,4)) NUS=3
          xgc=0.
          ygc=0.
          Do j=1,NUS
            xgc=xgc+x(KX(m,j))
            ygc=ygc+y(KX(m,j))
          Enddo
          TotalRootDensity=RUTDEN(m)+YRL(m)/Area(m)
          Write(82,140)t, Date,  xgc/NUS,ygc/NUS, m, RTWT(m),
     !     TotalRootDensity, SinkArr(m),cSinkArr(m),gSinkArr(m),Area(m)
        Enddo
C      Endif     
130   format('     Time           x             y     Temp   RTWT_',
     !        '          Sink_       cSink       gSink'/)
 140   format(1x,f14.6,A12,2(f11.2),i8,7(G11.3))
      Return
      End      
C/***************************************************/
      subroutine ArrNodOut(Date,x,y,NumNP,NumEl,t,NumElD,Head,      
     !          ThetaArr, vzArr,vxArr,QArr,ConcArr,TemprArr,gArr,Capa,
     !           Conduct)
      Character*10 Date
      Dimension Head(NumNP),x(NumNP),y(NumNP),
     !         ThetaArr(NumNP), vzArr(NumNP),vxArr(NumNP),QArr(NumNP),
     !         ConcArr(NumNP),TemprArr(NumNP),gArr(NumNP),Capa(NumNP),
     !         Conduct(NumNP)

        !MyDate=ToCalendarDate(mm,id,iyyy)
        Do m=1,NumNP
          Write(81,140)t, date, x(m),y(m),Head(m),ThetaArr(m),
     !     vzArr(m),vxArr(m),QArr(m),ConcArr(m),TemprArr(m),
     !      gArr(NumNP),Capa(m),Conduct(m)
        Enddo
130   format(//' Time            x            y       h____',
     !         '   Theta      vz____    vx____   Q____    Conc_ ',
     !         '    Tempr_     g____  Capa____')
 140   format(1x, f14.6,a12,2(f11.2),10(G11.3))
      Return
      End
  
  
