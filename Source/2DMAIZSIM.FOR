*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     2DSOIL  - Numerical model of two-dimensional soil environment    *
*                           for field crops                            *
*                                                                      *
*                              Release 1.1.4                              *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      Include 'public.ins'
      include 'puplant.ins'
	Include 'Puweath.ins'
	Include 'Nitvar.ins'
	
	character*132 RunDat
	integer len, status
	RunFile="Run.Dat"
	Call get_command_argument(1, RunDat, len, status)
	if (len.gt.0) then
	  RunFile=RunDat
	  endif
	  	
c	This is for testing
      Call Initialize()
      Call Get_Grid_and_Boundary()
1     Call Synchronizer()
      Call Crop()
C       
C SetSurfaceD --  reading daily weather
C SetSurfaceH --  reading hourly weather
      if ((HourlyWeather+DailyWeather.gt.1).or.
     !   ((HourlyWeather+DailyWeather).eq.0))then
              Write (*,*) "error in weather file type"
         
        Endif
      if (HourlyWeather.eq.1) Call SetSurfaceH()
       if (DailyWeather.eq.1) Call SetSurfaceD()                         
      Call SetTDB()
C      Call AutoIrrigate()
C      Call Mngm()

      Starter=1.0D0
C      Call CROP(PCRL,MSW1, NumNP,NumSol,Sink(1),
c     !      NumBP,TNext(1),NumMod,Starter)
C     Call Carbon_Partitioning()   
C      Call Rootgrow()
C      Call WaterUptake()  
C      Call SoluteUptake()
C      Call GasUptake()
      Call SetBoundary()
      Call Surface_Water_Balance_Adjustment()
      Call WaterMover()
      Call SoluteMover()
      Call HeatMover()
C      Call GasMover()
C      Call SoilNitrogen()
C      Call MacroChem()
      Call Water_Mass_Balance 
      if (outputSoilYes.gt.0) then 
         Call Output()
        endif

      
      Goto 1
      End         
