c ************************************************************************************************************
c There are three parts in the program, representing three trials we made for the runoff simulation
c Part I: Assume all of the runoff move to the next computation cell instantaneously. The Mass balance will work.
c Part II: Solve the fully coupled Saint-Venant, with some flat slope, the "h" will increase to a unreasonable value, so the diff
c         equation system may not always work
c Part III: Not only consider the coupled Saint-Venant, but also set a threshold for "h", the rest will go to the runoff.
c         The momentun equation will be used to estimate the flux speed.???
c ************************************************************************************************************

      Subroutine Surface_Water_Balance_Adjustment()
      include 'Public.ins'
      include 'puplant.ins'
      
      
      double precision Slope_0_R, n_Stiff, g_accel, Slope_f
c SURNODE to index the node order from the whole node array      
c SURNODE_Sur to index the node order from the boundary node array  
      integer  FluxDir(15),iteration_Flux, iteration_Head
      double precision RunoffLeft_temp, RunoffRight_temp, iteration_Num
      Common /SurWaterBalance/ModNum,h_Stay(NumNPD),q_Flux(NumNPD), 
     !       h_Stay_temp(NumNPD),q_Flux_temp(NumNPD),R_I(NumNPD),
     !       h_Stay_Old(NumNPD),q_Flux_Old(NumNPD),CriticalH_R,
     !       SlopeHV(NumNPD,2),q_Flux_Node(NumNPD),
     !       time_sur,h_Stay_ini,RunoffLeft,RunoffRight,
     !       RunoffLeft_old, RunoffRight_old,
     !       FluxLimit, HeadLimit
      
      If (lInput.eq.1) then
c ****************** This part is an extension of hour.wea ******************
          CriticalH_R=1.0D0
          h_Stay_ini=0.0D0
          
          ModNum=1
          Update=0
          do i=1,NumNPD
              RO(i)=0.0D0
              VarBW_Old_1(i)=0.0D0
              VarBW_Old_2(i)=0.0D0
              VarBW_Old_3(i)=0.0D0
          enddo
          
          do i=1,NumNPD
              h_Stay(i)=h_Stay_ini
              q_Flux(i)=0.0D0
              h_Stay_temp(i)=h_Stay_ini
              q_Flux_temp(i)=0.0D0
              q_Flux_Node(i)=0.0D0
              R_I(i)=0.0D0
          enddo
          
c Part IV. Use the New Grid. We use a step-by-step method to read the grid on the surface.
c          As well as determine the horizontal distance, vertical distance between each node.
c ******************************************************************************************  
          SurIndex=0
          do i=1,NumBP
              n=KXB(i)
              k=CodeW(n)
              If(K.eq.4.or.K.eq.-4) then
                  SurIndex=SurIndex+1
                  SURNODE(SurIndex)=n
              endif
          enddo
          
          InElement=0          
          do i=1,SurIndex
             do j=1,NumEL
               do kk=1,4
                  if(SURNODE(i).eq.KX(j,kk)) then
                      InElement=1
                  endif
                enddo
                if(InElement.eq.1) then
                  do kk=1,4
                  if(CodeW(KX(j,kk)).eq.4.or.CodeW(KX(j,kk)).eq.-4) then
                     if(i.eq.1) then
                        if(KX(j,kk).ne.SURNODE(i)) then
                          SURNODE(i+1)=KX(j,kk)
                        endif
                     else
                        if(KX(j,kk).ne.SURNODE(i).
     &                    and.KX(j,kk).ne.SURNODE(i-1)) then
                           SURNODE(i+1)=KX(j,kk) 
                        endif
                     endif
                    endif
                  enddo
                 endif
                 InElement=0
               enddo
          enddo
          
          do i=1,SurIndex
              SlopeHV(i,1)=X(SURNODE(i))
              SlopeHV(i,2)=Y(SURNODE(i))
          enddo
          
          do n=1,SurIndex
              do i=1,NumBP
                  if(SURNODE(n).eq.KXB(i)) then
                        SURNODE_Sur(n)=i
                  endif
              enddo
          enddo
          
          RunoffLeft=0.0D0
          RunoffRight=0.0D0 
          RunoffLeft_temp=0.0D0
          RunoffRight_temp=0.0D0
          RunoffLeft_old=0.0D0
          RunoffRight_old=0.0D0
c ************************************
c In the iteration we will not care about the error if the magnitude of the flux or head is below the limiter
          FluxLimit=0.00000000001D0
          HeadLimit=0.00000000001D0
          
      Endif
            
      If(Update.ge.1.or.lInput.eq.1) then
          do i=1,NumBP
              VarBW_Old_1(i)=VarBW(i,1)
              VarBW_Old_2(i)=VarBW(i,2)
              VarBW_Old_3(i)=VarBW(i,3)
          enddo
c ************************************ 
c We can play with the rainfall rate here. But we force the surface nodes from 2 to 7
c Manually change the rainfall, can be eliminated with real rainfall condition
c ************************************
c          do i=2,7
c              VarBW_Old_1(i)=VarBW(i,1)*0.40D0
c          enddo
c ************************************************************************
          Update=0
      endif
      
c Part I. Just move the runoff to the next cell.      
c ****************** Simple equation for moving runoff ******************  
c          do i=1,NumBP
c              n=KXB(i)
c              k=CodeW(n)
c              If(K.eq.4.or.K.eq.-4) then 
c                  If (i.gt.1) then
c                      VarBW(i,1)=VarBW_Old_1(i)
c     &                      +max(RO(i-1)/Width(i),0.0d0)
c                      VarBW(i,3)=VarBW_Old_2(i)-VarBW(i,1)
c                      Q(n)=-Width(i)*VarBW(i,3)
c                  
c                      If (Q(n).gt.0.0) then
c                          CodeW(n)=-4
c                      endif
c                  endif
c              endif
c          enddo
c ******************************************************************************************
C Part I. End

c Part II. Solve the fully coupled PDE
c Problems: the "h" values will increase up high, and it cannot handle the decreasing of rainfall
c ****************** Begin to write a new code based on Saint-Venant Eqs, refer to Kouznetsov et al., 2007 ******************
c          do i=1,NumBP
c              n=KXB(i)
c              k=CodeW(n)
c              If(K.eq.4.or.K.eq.-4) then 
c                  h_Stay_Old(i)=h_Stay(i)
c                  q_Flux_Old(i)=q_Flux(i)
c ****************** use iteration to partition the real runoff and amount of water stay
c                  R_I(i)=RO(i)/width(i)
c                  n_Stiff=0.0D0
c                  Slope_0=0.01D0
c                  g_accel=0.0113425926D0
c ****************** We need some hard judgement for the unit, we use cm and day
c ****************** Explicit Version
c10                Slope_f=((n_Stiff*q_Flux(i))**2.0D0)
c     &                /(h_Stay(i)**(3.333D0))
c                  if (i.eq.1) then
c                      h_Stay_temp(i)=h_Stay_Old(i)
c     &                   -step*q_Flux(i)/width(i)+R_I(i)*step
c                      q_Flux_temp(i)=q_Flux_Old(i)
c     &                   -step*(q_Flux(i)*q_Flux(i)/h_Stay(i)
c     &                   +0.5D0*g_accel*h_Stay(i)*h_Stay(i))/width(i)
c     &                   +g_accel*h_Stay(i)*(Slope_0-Slope_f)*step
c                  else
c                       h_Stay_temp(i)=h_Stay_Old(i)
c     &                      -step*q_Flux(i)/width(i)
c     &                      +R_I(i)*step
c                       q_Flux_temp(i)=q_Flux_Old(i)
c     &                   -step*(q_Flux(i)*q_Flux(i)/h_Stay(i)
c     &                   +0.5D0*g_accel*h_Stay(i)*h_Stay(i))/width(i)
c     &                   +step*(q_Flux(i-1)*q_Flux(i-1)/h_Stay(i-1)
c     &                  +0.5D0*g_accel*h_Stay(i-1)*h_Stay(i-1))/width(i)
c     &                   +g_accel*h_Stay(i)*(Slope_0-Slope_f)*step
c                      
c                  endif
c                  q_Flux_temp(i)=max(q_Flux_temp(i),0.0D0)
c                  h_Stay_temp(i)=max(h_Stay_temp(i),h_Stay_ini)
c                  if (h_Stay_temp(i).le.h_Stay_ini) then
c                      q_Flux_temp(i)=0.0D0                      
c                  endif
c                                    
c                  if (abs(q_Flux_temp(i)-q_Flux(i)).le.0.0001*q_Flux(i)
c     &             .and.abs(h_Stay_temp(i)-h_Stay(i))
c     &             .le.0.0001*h_Stay(i)) then
c                          q_Flux(i)=q_Flux_temp(i)
c                          h_Stay(i)=h_Stay_temp(i)
c                  else
c                          q_Flux(i)=q_Flux_temp(i)
c                          h_Stay(i)=h_Stay_temp(i)
c                          goto 10
c                  endif
c                                               
c                  if(i.eq.1) then
c                      VarBW(i,1)=VarBW_Old_1(i)
c                  else
c                      VarBW(i,1)=VarBW_Old_1(i)+q_Flux(i-1)/width(i)
c                  endif
c                  
c                  VarBW(i,3)=VarBW_Old_2(i)-VarBW(i,1)
c                  Q(n)=-Width(i)*VarBW(i,3)
c                                               
c                  If (Q(n).gt.0.0) then
c                      CodeW(n)=-4
c                  endif
c                  
c ****************** Temperatory Check ******************                
c                  h_stay_water=h_stay_water
c     &                 +(h_Stay(i)-h_Stay_Old(i))*width(i)
c ******************************************************
c
c              endif
c          enddo
c ******************************************************************************************
c Part II. End


c Part III. Add a threshold for "h", then solve the fully coupled PDE
c ******************************************************************************************
c          do i=1,NumBP
c              n=KXB(i)
c              k=CodeW(n)
c              If(K.eq.4.or.K.eq.-4) then 
c                  h_Stay_Old(i)=h_Stay(i)
c                  q_Flux_Old(i)=q_Flux(i)
c ****************** use iteration to partition the real runoff and amount of water stay
c                  R_I(i)=RO(i)/width(i)
c                  n_Stiff=0.0D0
c                  Slope_0=0.01D0
c                  g_accel=0.0113425926D0
c ****************** We need some hard judgement for the unit, we use cm and day
c ****************** Explicit Version
c10                Slope_f=((n_Stiff*q_Flux(i))**2.0D0)
c     &                /(h_Stay(i)**(3.333D0))
c                  if (i.eq.1) then
c                      
c                      q_Flux_temp(i)=q_Flux_Old(i)
c     &                   -step*(q_Flux(i)*q_Flux(i)/h_Stay(i)
c     &                   +0.5D0*g_accel*h_Stay(i)*h_Stay(i))/width(i)
c     &                   +g_accel*h_Stay(i)*(Slope_0-Slope_f)*step
c                      h_Stay_temp(i)=h_Stay_Old(i)
c     &                   -step*q_Flux(i)/width(i)+R_I(i)*step
c     
c                  else
c                       
c                      q_Flux_temp(i)=q_Flux_Old(i)
c     &                   -step*(q_Flux(i)*q_Flux(i)/h_Stay(i)
c     &                   +0.5D0*g_accel*h_Stay(i)*h_Stay(i))/width(i)
c     &                   +step*(q_Flux(i-1)*q_Flux(i-1)/h_Stay(i-1)
c     &                  +0.5D0*g_accel*h_Stay(i-1)*h_Stay(i-1))/width(i)
c     &                   +g_accel*h_Stay(i)*(Slope_0-Slope_f)*step
c                      h_Stay_temp(i)=h_Stay_Old(i)
c     &                   -step*q_Flux(i)/width(i)
c     &                   +R_I(i)*step
c                      
c                  endif
c                  
c                  h_Stay_temp(i)=min(max(h_Stay_temp(i),
c     &                 h_Stay_ini),CriticalH_R)
c ************************************************************************
c If "h" is smaller than CriticalH_R, the following updates will not change the flux value
c If "h" is larger than CriticalH_R, the following updates will change the flux value using the conservation of mass
c Will the update flux value work for the conservation of momentun?
c ************************************************************************
c                  q_Flux_temp(i)=R_I(i)*width(i)
c     &                 +(h_Stay_Old(i)-h_Stay_temp(i))*width(i)/step
c                  q_Flux_temp(i)=max(q_Flux_temp(i),0.0D0)
c                  
c                  if (h_Stay_temp(i).le.h_Stay_ini) then
c                      q_Flux_temp(i)=0.0D0                      
c                  endif
c                                    
c                  if (abs(q_Flux_temp(i)-q_Flux(i)).le.0.0001*q_Flux(i)
c     &             .and.abs(h_Stay_temp(i)-h_Stay(i))
c     &             .le.0.0001*h_Stay(i)) then
c                         q_Flux(i)=q_Flux_temp(i)
c                          h_Stay(i)=h_Stay_temp(i)
c                  else
c                          q_Flux(i)=q_Flux_temp(i)
c                         h_Stay(i)=h_Stay_temp(i)
c                          goto 10
c                  endif
c                                               
c                  if(i.eq.1) then
c                      VarBW(i,1)=VarBW_Old_1(i)
c                  else
c                      VarBW(i,1)=VarBW_Old_1(i)+q_Flux(i-1)/width(i)
c                  endif
c                  
c                  VarBW(i,3)=VarBW_Old_2(i)-VarBW(i,1)
c                  Q(n)=-Width(i)*VarBW(i,3)
c                                               
c                  If (Q(n).gt.0.0) then
c                      CodeW(n)=-4
c                  endif
c                  
c ****************** Temperatory Check                
c                  h_stay_water=h_stay_water
c     &                 +(h_Stay(i)-h_Stay_Old(i))*width(i)
c ******************
c
c              endif
c          enddo
c ******************************************************************************************
c Part III. End


c Part IV. Use the New Grid. Add a threshold for "h", then solve the fully coupled PDE
c ******************************************************************************************          
          

          
          do n=1,SurIndex
              i=SURNODE_Sur(n)
              k=SURNODE(n)
              h_Stay_Old(n)=h_Stay(n)
              q_Flux_Old(n)=q_Flux(n)
              R_I(n)=max(RO(k)/width(i),0.0D0)
          enddo
          RunoffLeft_old=RunoffLeft
          RunoffRight_old=RunoffRight
          
c ************ Set the boundary H thresholds at 1.0D0
          iteration_Num=0          

104       iteration_Flux=0
          iteration_Head=0
          
          RunoffLeft_temp=0.0D0
          RunoffRight_temp=0.0D0
         
          n_Stiff=1.0D0
          g_accel=0.0113425926D0
c ****************** We need some hard judgement for the unit, we use cm and day  

          do i=1,NumNPD
              h_Stay_temp(i)=h_Stay_ini
              q_Flux_temp(i)=0.0D0
              q_Flux_Node(i)=0.0D0
          enddo
          
          do n=1,SurIndex-1
c ****************** use iteration to partition the real runoff and amount of water stay 
c ****************** take the absolute value of the slope to ensure the following calculation
              Slope_0_R=abs((SlopeHV(n+1,2)-SlopeHV(n,2))
     &                   /(SlopeHV(n,1)-SlopeHV(n+1,1)))

c ****************** Explicit Version
              if(q_Flux(n+1).eq.0) then
                  Slope_f_1=0.0D0
              else
                  Slope_f_1=((n_Stiff*q_Flux(n))**2.0D0)
     &             /(h_Stay(n+1)**(3.333D0))
              endif
              if(q_Flux(n).eq.0) then
                  Slope_f_2=0.0D0
              else
                  Slope_f_2=((n_Stiff*q_Flux(n)**2.0D0)
     &            /(h_Stay(n)**(3.333D0)))
              endif
              Slope_f=0.50D0*(Slope_f_1+Slope_f_2)
              
            if(n.eq.1) then
               if((SlopeHV(n,2)+h_Stay(n)).ge.
     &            (SlopeHV(n+1,2)+h_Stay(n+1))) then
c Water will flow from 1 to 2
c Always follow the upwind direction                      
                    q_Flux_temp(n)=q_Flux_Old(n)
     &                -step*q_Flux(n)*q_Flux(n)/h_Stay(n)
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                +g_accel*h_Stay(n)
     &                *max(Slope_0_R-Slope_f,0.0D0)*step
                    q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
                    
                    RunoffLeft_temp=0.0D0
               else
c Water will flow from 2 to 1, momentun source of 2 from 3?
c Always follow the upwind direction                      
                  if ((SlopeHV(n+1,2)+h_Stay(n+1)).ge.
     &                   (SlopeHV(n+2,2)+h_Stay(n+2))) then          
                    q_Flux_temp(n)=q_Flux_Old(n)
     &                   +step*(q_Flux(n)*q_Flux(n)/h_Stay(n+1))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                -g_accel*h_Stay(n+1)
     &                *max(Slope_0_R-Slope_f,0.0D0)*step
                     q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
                  else
                     q_Flux_temp(n)=q_Flux_Old(n)
     &                   +step*(q_Flux(n)*q_Flux(n)/h_Stay(n+1)
     &                     -q_Flux(n+1)*q_Flux(n+1)/h_Stay(n+2))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                -g_accel*h_Stay(n+1)
     &                *max(Slope_0_R-Slope_f,0.0D0)*step 
                     q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
                  endif
c first calculate the runoff from the left edge
                  RunoffLeft_temp=RunoffLeft_old
     &               +step*(RunoffLeft*RunoffLeft/h_Stay(n)
     &                     -q_Flux(n)*q_Flux(n)/h_Stay(n+1))
     &                     /width(SURNODE_sur(n))
     &                  +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                  /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   -g_accel*h_Stay(n)
     &                 *max(Slope_0_R-Slope_f,0.0D0)*step 
                  RunoffLeft_temp=min(RunoffLeft_temp,0.0D0)
               endif
               
             elseif(n.gt.1.and.n.lt.SurIndex-1) then
               if((SlopeHV(n,2)+h_Stay(n)).ge.
     &            (SlopeHV(n+1,2)+h_Stay(n+1))) then
c Water will flow from n to n+1, momentun source of n-1 from n?
c Always follow the upwind direction                  
                  if ((SlopeHV(n,2)+h_Stay(n)).ge.
     &                   (SlopeHV(n-1,2)+h_Stay(n-1))) then          
                     q_Flux_temp(n)=q_Flux_Old(n)
     &                   -step*(q_Flux(n)*q_Flux(n)/h_Stay(n))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                +g_accel*h_Stay(n)
     &                 *max(Slope_0_R-Slope_f,0.0D0)*step
                     q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
                  else
                     q_Flux_temp(n)=q_Flux_Old(n)
     &                   -step*(q_Flux(n)*q_Flux(n)/h_Stay(n)
     &                     -q_Flux(n-1)*q_Flux(n-1)/h_Stay(n-1))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                +g_accel*h_Stay(n)
     &                 *max(Slope_0_R-Slope_f,0.0D0)*step  
                      q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
                  endif
               else
c Water will flow from n+1 to n, momentun source of n+2 from n+1?
c Always follow the upwind direction                     
                  if ((SlopeHV(n+1,2)+h_Stay(n+1)).ge.
     &                   (SlopeHV(n+2,2)+h_Stay(n+2))) then          
                     q_Flux_temp(n)=q_Flux_Old(n)
     &                   +step*(q_Flux(n)*q_Flux(n)/h_Stay(n+1))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                -g_accel*h_Stay(n+1)
     &                *max(Slope_0_R-Slope_f,0.0D0)*step
                     q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
                  else
                     q_Flux_temp(n)=q_Flux_Old(n)
     &                   +step*(q_Flux(n)*q_Flux(n)/h_Stay(n+1)
     &                     -q_Flux(n+1)*q_Flux(n+1)/h_Stay(n+2))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                -g_accel*h_Stay(n+1)
     &                *max(Slope_0_R-Slope_f,0.0D0)*step
                     q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
                  endif
               endif
             else
c Finish the calculation for n=SurIndex-1     
               if((SlopeHV(n,2)+h_Stay(n)).ge.
     &            (SlopeHV(n+1,2)+h_Stay(n+1))) then
c water flow to the end point, need to calculate the momentun from n-2
                  if ((SlopeHV(n,2)+h_Stay(n)).ge.
     &                   (SlopeHV(n-1,2)+h_Stay(n-1))) then   
                      q_Flux_temp(n)=q_Flux_Old(n)
     &                   -step*(q_Flux(n)*q_Flux(n)/h_Stay(n))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                +g_accel*h_Stay(n)
     &                *max(Slope_0_R-Slope_f,0.0D0)*step
                     q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
                  else
                     q_Flux_temp(n)=q_Flux_Old(n)
     &                   -step*(q_Flux(n)*q_Flux(n)/h_Stay(n)
     &                     -q_Flux(n-1)*q_Flux(n-1)/h_Stay(n-1))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                +g_accel*h_Stay(n)
     &                *max(Slope_0_R-Slope_f,0.0D0)*step
                     q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
                  endif
                  RunoffRight_temp=RunoffRight_old
     &               -step*(RunoffRight*RunoffRight/h_Stay(n+1)
     &                     -q_Flux(n)*q_Flux(n)/h_Stay(n+1))
     &                     /width(SURNODE_sur(n+1))
     &                  +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                  /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +g_accel*h_Stay(n+1)
     &                 *max(Slope_0_R-Slope_f,0.0D0)*step 
                  RunoffRight_temp=max(RunoffRight_temp,0.0D0)
                else
c water flow to from end point, to SurIndex-1
                   q_Flux_temp(n)=q_Flux_Old(n)
     &                   +step*(q_Flux(n)*q_Flux(n)/h_Stay(n+1))
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                   +0.5D0*g_accel*(h_Stay(n)*h_Stay(n)
     &                    -h_Stay(n+1)*h_Stay(n+1))*step 
     &                /abs(SlopeHV(n,1)-SlopeHV(n+1,1))
     &                -g_accel*h_Stay(n+1)
     &                *max(Slope_0_R-Slope_f,0.0D0)*step
                     q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
                     
                     RunoffRight_temp=0.0D0
                endif
             endif
             
        enddo
c finishe the flux estimation, now go the water head estimation
c The following steps will make the h_temp goes infinitely high

c        do n=1,SurIndex
c           i=SURNODE_Sur(n)
c           if(n.eq.1) then
c              h_Stay_temp(n)=h_Stay_Old(n)
c     &          -step*q_Flux(n)/width(i)+R_I(n)*step
c           elseif(n.gt.1.and.n.lt.SurIndex) then
c               h_Stay_temp(n)=h_Stay_Old(n)-step*q_Flux(n)/width(i)
c     &            +step*q_Flux(n-1)/width(i)+R_I(i)*step
c           else
c               h_Stay_temp(n)=h_Stay_Old(n)
c     &            +step*q_Flux(n-1)/width(i)+R_I(i)*step
c           endif
c           h_Stay_temp(n)=max(h_Stay_temp(n),
c     &                 h_Stay_ini)
c        enddo

c further adjust for the surface ponded height based on the threshold
c        print*,'f',n,'\',q_Flux_temp(n),'\',q_Flux(n)
        do n=1,SurIndex
           i=SURNODE_Sur(n)
           k=SURNODE(n)
           if(n.eq.1) then
              h_Stay_temp(n)=h_Stay_Old(n)
     &          -max(step*q_Flux_temp(n)/width(i),0.0D0)
     &          +min(step*RunoffLeft_temp/width(i),0.0D0)+R_I(n)*step
              h_Stay_temp(n)=max(h_Stay_temp(n),h_Stay_ini)
              
               if((SlopeHV(n,2)+h_Stay(n)).le.
     &            (SlopeHV(n+1,2)+h_Stay(n+1))) then
                    if(h_Stay_temp(n).gt.CriticalH_R) then
c                      h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
c                      RunoffLeft_temp=-(R_I(n)*width(i)
c     &                 +(h_Stay_Old(n)-h_Stay_temp(n))*width(i)/step)
                       RunoffLeft_temp=RunoffLeft_temp
     &                    -(h_Stay_temp(n)-CriticalH_R)*width(i)/step
                       h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
                     endif
                     
                else
                    if(h_Stay_temp(n).gt.CriticalH_R) then
c                      h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
c                     q_Flux_temp(n)=R_I(n)*width(i)
c     &                 +(h_Stay_Old(n)-h_Stay_temp(n))*width(i)/step

                      q_Flux_temp(n)=q_Flux_temp(n)
     &                  +(h_Stay_temp(n)-CriticalH_R)*width(i)/step
                      h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
                       
                      RunoffLeft_temp=0.0D0
                     endif           
c                   q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
                endif
              
            elseif(n.gt.1.and.n.lt.SurIndex) then
               h_Stay_temp(n)=h_Stay_Old(n)
     &            -max(step*q_Flux_temp(n)/width(i), 0.0D0)
     &            +min(step*q_Flux_temp(n-1)/width(i), 0.0D0)
     &            +R_I(n)*step
               h_Stay_temp(n)=max(h_Stay_temp(n),h_Stay_ini)
               if((SlopeHV(n,2)+h_Stay(n)).ge.
     &            (SlopeHV(n-1,2)+h_Stay(n-1))) then
c               if(q_Flux(n-1).le.0) then
c Previuos time step will not take care of q_Flux(n-1)
                   if((SlopeHV(n,2)+h_Stay(n)).le.
     &                  (SlopeHV(n+1,2)+h_Stay(n+1))) then
c                   if(q_Flux(n).le.0) then
                     if(h_Stay_temp(n).gt.CriticalH_R) then
c                       h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
c                       q_Flux_temp(n-1)=-R_I(n)*width(i)
c     &                    -(h_Stay_Old(n)-h_Stay_temp(n))*width(i)/step
                        q_Flux_temp(n-1)=q_Flux_temp(n-1)
     &                     -(h_Stay_temp(n)-CriticalH_R)*width(i)/step
                        h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
                        
                     endif
                   else
                     if(h_Stay_temp(n).gt.CriticalH_R) then
c                       h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
c                       q_Flux_temp(n-1)=R_I(n)*width(i)
c     &                   +(h_Stay_Old(n)-h_Stay_temp(n))*width(i)/step
c redistribute the flux based on the slope
c                       Slope_0_n_1=abs((SlopeHV(n-1,2)-SlopeHV(n,2))
c     &                   /(SlopeHV(n,1)-SlopeHV(n-1,1)))
c                       Slope_0_n=abs((SlopeHV(n+1,2)-SlopeHV(n,2))
c     &                   /(SlopeHV(n,1)-SlopeHV(n+1,1)))
c                       Fraction_n_1=Slope_0_n_1/(Slope_0_n_1+Slope_0_n)
c                       Fraction_n=Slope_0_n/(Slope_0_n_1+Slope_0_n)
c                       q_Flux_temp(n-1)=-Fraction_n_1*q_Flux_temp(n-1)
c                       q_Flux_temp(n)=Fraction_n*q_Flux_temp(n-1)
                       Add_Flux=(h_Stay_temp(n)-CriticalH_R)
     &                        *width(i)/step
                       Slope_0_n_1=abs((SlopeHV(n-1,2)-SlopeHV(n,2))
     &                     /(SlopeHV(n,1)-SlopeHV(n-1,1)))
                       Slope_0_n=abs((SlopeHV(n+1,2)-SlopeHV(n,2))
     &                     /(SlopeHV(n,1)-SlopeHV(n+1,1)))
                       Fraction_n_1=Slope_0_n_1/(Slope_0_n_1+Slope_0_n)
                       Fraction_n=Slope_0_n/(Slope_0_n_1+Slope_0_n)
                       
                       q_Flux_temp(n-1)=q_Flux_temp(n-1)
     &                        -Fraction_n_1*Add_Flux
                       q_Flux_temp(n)=q_Flux_temp(n)
     &                        +Fraction_n*Add_Flux
                       h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
                    endif
                  endif
                else
c The q_Flux(n-1) are taken cared by previous steps, then just calculate new q_Flux(n)
c if and only if q_Flux(n)>0
                   if((SlopeHV(n,2)+h_Stay(n)).ge.
     &                  (SlopeHV(n+1,2)+h_Stay(n+1))) then
c                   if(q_Flux(n).ge.0) then
                     if(h_Stay_temp(n).gt.CriticalH_R) then
c                       h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
c                       q_Flux_temp(n)=R_I(n)*width(i)
c     &                    +(h_Stay_Old(n)-h_Stay_temp(n))*width(i)/step
                       
                       q_Flux_temp(n)=q_Flux_temp(n)
     &                    +(h_Stay_temp(n)-CriticalH_R)*width(i)/step
                       h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
                       
                     endif
                   else
                     if(h_Stay_temp(n).gt.CriticalH_R) then
                       
                       q_Flux_Node(n)=(h_Stay_temp(n)-CriticalH_R)
     &                     *width(i)/step
                       h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
                       
                     endif
                   endif
c if q_Flux(n).le.0, it will be calculate in the next loop
                endif
                
            else
            
                h_Stay_temp(n)=h_Stay_Old(n)
     &            -max(step*RunoffRight_temp/width(i),0.0D0)
     &            +min(step*q_Flux_temp(n-1)/width(i),0.0D0)
     &            +R_I(n)*step
                h_Stay_temp(n)=max(h_Stay_temp(n),h_Stay_ini)

                if((SlopeHV(n,2)+h_Stay(n)).le.
     &             (SlopeHV(n-1,2)+h_Stay(n-1))) then
c                  if(q_Flux(n-1).ge.0) then
                     if(h_Stay_temp(n).gt.CriticalH_R) then
c                       h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
c                       RunoffRight=R_I(n)*width(i)
c     &                  +(h_Stay_Old(n)-h_Stay_temp(n))*width(i)/step
                        RunoffRight_temp=RunoffRight_temp
     &                   +(h_Stay_temp(n)-CriticalH_R)*width(i)/step
                        h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
                     endif
                else
                     if(h_Stay_temp(n).gt.CriticalH_R) then
c                        h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
c                        q_Flux_temp(n-1)=-R_I(n)*width(i)
c     &                    -(h_Stay_Old(n)-h_Stay_temp(n))*width(i)/step\
                        
                        q_Flux_temp(n-1)=q_Flux_temp(n-1)
     &                    -(h_Stay_temp(n)-CriticalH_R)*width(i)/step
                        h_Stay_temp(n)=min(h_Stay_temp(n),CriticalH_R)
                        
                        RunoffRight_temp=0.0D0
                     endif
                endif
            endif
         enddo
c finish update of h, need further work for setting the critical H value for each node.

c calculate the runoff at two edge
c left edge for n=1

                  
c ************************************************************************
c If "h" is smaller than CriticalH_R, the following updates will not change the flux value
c If "h" is larger than CriticalH_R, the following updates will change the flux value using the conservation of mass
c Will the update flux value work for the conservation of momentun?
c ************************************************************************
c                  q_Flux_temp(i)=R_I(i)*width(i)
c     &                 +(h_Stay_Old(i)-h_Stay_temp(i))*width(i)/step
c                  q_Flux_temp(i)=max(q_Flux_temp(i),0.0D0)
                  
c                  if (h_Stay_temp(i).le.h_Stay_ini) then
c                      q_Flux_temp(i)=0.0D0                      
c                  endif
c Error esitimation for iteration
          
          iteration_Flux=0
          iteration_Head=0
          iteration_Num=iteration_Num+1
c          iteration_Num=1.0D0
          do n=1,SurIndex
              q_Flux_temp(n)=q_Flux_temp(n)/iteration_Num+
     &           (iteration_Num-1)*q_Flux(n)/iteration_Num
              h_Stay_temp(n)=h_Stay_temp(n)/iteration_Num+
     &           (iteration_Num-1)*h_Stay(n)/iteration_Num
          enddo
          RunoffRight_temp=RunoffRight_temp/iteration_Num+
     &          (iteration_Num-1)*RunoffRight/iteration_Num
          RunoffLeft_temp=RunoffLeft_temp/iteration_Num+
     &          (iteration_Num-1)*RunoffLeft/iteration_Num
          
          
          do n=1,SurIndex
             if(abs(q_Flux_temp(n)-q_Flux(n)).gt.
     &             0.01*abs(q_Flux(n))) then
               if(abs(q_Flux(n)).gt.FluxLimit) then
                  iteration_Flux=1
                endif
c               print*,'f',n,'\',q_Flux_temp(n),'\',q_Flux(n)
             endif
                          
             if(abs(h_Stay_temp(n)-h_Stay(n)).gt.
     &             0.01*abs(h_Stay(n))) then
               if(abs(h_Stay(n)).gt.HeadLimit) then
                  iteration_Head=1
                endif
c               print*,'h',n,'\',h_Stay_temp(n)-h_Stay(n)
             endif
          enddo
          
          if(abs(RunoffLeft_temp-RunoffLeft).gt.
     &               0.01*abs(RunoffLeft)) then
               if(abs(RunoffLeft).gt.FluxLimit) then
                  iteration_Flux=1
                endif
          endif
             
          if(abs(RunoffRight_temp-RunoffRight).gt.
     &               0.01*abs(RunoffRight)) then
               if(abs(RunoffRight).gt.FluxLimit) then
                  iteration_Flux=1
                endif
          endif

          if (iteration_Flux.eq.1.or.iteration_Head.eq.1) then
c            Iteration_Num=Iteration_Num+1
            do n=1,SurIndex
              q_Flux(n)=q_Flux_temp(n)
              h_Stay(n)=h_Stay_temp(n)
            enddo
            RunoffLeft=RunoffLeft_temp
            RunoffRight=RunoffRight_temp
            goto 104
          else
            do n=1,SurIndex
              q_Flux(n)=q_Flux_temp(n)
              h_Stay(n)=h_Stay_temp(n)
            enddo
            RunoffLeft=RunoffLeft_temp
            RunoffRight=RunoffRight_temp
          endif

c Now assign the flux as water input at each node
          do n=1,SurIndex
             if (n.eq.1) then
                q_Flux_Node(n)=q_Flux_Node(n)-min(q_Flux(n),0.0D0)
     &               +max(RunoffLeft,0.0D0)
             elseif (n.gt.1.and.n.lt.SurIndex) then
                q_Flux_Node(n)=q_Flux_Node(n)+max(q_Flux(n-1),0.0D0)
     &               -min(q_Flux(n),0.0D0)
             else
                q_Flux_Node(n)=q_Flux_Node(n)+max(q_Flux(n-1),0.0D0)
     &               -min(RunoffRight,0.0D0)
             endif
          enddo
          
          do n=1,SurIndex
              k=SURNODE(n)            ! for index in the whole node set
              i=SURNODE_Sur(n)        ! for index in the boundary node set
              VarBW(i,1)=VarBW_Old_1(i)+q_Flux_Node(n)/width(i)
              VarBW(i,3)=VarBW_Old_2(i)-VarBW(i,1)
              Q(k)=-Width(i)*VarBW(i,3)
              If (Q(k).gt.0.0) then
                      CodeW(k)=-4
              endif
          enddo
                                                                 
c ******************************************************************************************
c Part IV. End

c          zhuangji=1.0D0
         print*,time,'\',h_Stay(3),q_Flux(3),RunoffLeft, RunoffRight
         Runoff_Terminal_Flux=Runoff_Terminal_Flux
     &                +(RunoffLeft+RunoffRight)*step
     
       if(time.ge.38751.0) then
           print*, 'ZZZ'
       endif
       
      return
      
      end