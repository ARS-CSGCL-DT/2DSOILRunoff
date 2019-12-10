# 2DSOILRunoff
Soil water simulation including infiltration and runoff based on 2DSOIL
2DSOILRunoff is a finite element model of soil water dynamics that was originally developed from SWMS_2D (Jerka Simunek and Tomas Vogel). This adaptation has been coded by Zhuangji Wang, Mikhail Kouznetsov and Dennis Timlin. It includes a runoff model based on the St. Venant's equation and uses the Heaviside step function to switch boundary conditions. It was branched from MAIZSIM, a mechanistic simulation model of maize growth and development. The version in this repository is a stand alone version without a crop model. We plan to reincorporated the plant model later in the future thus some elements of the plant model are still in the code. This mainly includes reading plant related files in the main input file.