c--- Grid of z-vaules for fragmentation functions

      integer num_z
      parameter(num_z = 45)
      
      double precision z_grid(num_z) 
      
      data z_grid  / 
     . +0.4627D-01,0.5131D-01,0.5690D-01,0.6311D-01,0.6999D-01,           
     . +0.7762D-01,
     . +0.8608D-01,0.9547D-01,0.1059D+00,0.1174D+00,0.1302D+00,
     . +0.1444D+00,
     . +0.1602D+00,0.1776D+00,0.1970D+00,0.2185D+00,0.2423D+00,
     . +0.2687D+00,
     . +0.2980D+00,0.3305D+00,0.3666D+00,0.4065D+00,0.4508D+00,
     . +0.5000D+00,
     . +0.5321D+00,0.5643D+00,0.5964D+00,0.6286D+00,0.6607D+00,
     . +0.6929D+00,
     . +0.7250D+00,0.7571D+00,0.7893D+00,0.8214D+00,0.8536D+00,
     . +0.8857D+00,
     . +0.9000D+00,0.9179D+00,0.9300D+00,0.9400D+00,0.9500D+00,
     . +0.9600D+00,      
     . +0.9700D+00,0.9800D+00,0.99D+00 /

c--- Grid of M**2 Values for fragmentation functions

      integer num_M2_4Flav, num_M2_5Flav
      parameter (num_M2_4Flav=8,num_M2_5Flav=22)
      
      double precision M2_4Flav(num_M2_4Flav), M2_5Flav(num_M2_5Flav)
      
      data M2_4Flav / 
     . 2.D0,2.8D0,3.9D0,5.4D0,7.5D0,10.5D0,14.5D0,20.24D0/
      
      data M2_5Flav /
     . +20.26D0,33.D0,55.D0,93.D0,156.D0,261.D0,
     . +435.D0,732.D0,1224.D0,2048.D0,
     . +3427.D0,5733.D0,9592.D0,16047.D0,
     .  +26847.D0,44915.D0,75144.D0,1.257d+5,
     . +2.1d+5,3.52d+5,5.89d+5,9.848d+5 /

      


      
     
