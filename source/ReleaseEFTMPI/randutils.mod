  �  8   k820309    `          17.0        ��[                                                                                                           
       RandUtils.f90 RANDUTILS                                                     
                                                              u #RANDROTATIONS    #RANDROTATIOND                                                      KIND #         @                                                      #MSG              
                                                    1 #         @                                                       #THIS                                                                  #TTIMER 	   %         @                               
                    
       #THIS                                                                  #TTIMER 	   #         @                                                       #THIS    #MSG    #START                                                                  #TTIMER 	             
                                                    1                                                
                                                                                                                                                                                                                             #         @      X                                                 #R    #N             D                                                     	       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                
                                             #         @      X                                                 #R    #N             D                                                     
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                
                                             #         @                                                       #I    #I2              
 @                                                    
 @                                          #         @                                                     #RANDUTILS!RMARIN!RASET1    #IJ "   #KL #                                                  (                    #RMARIN%RASET1%U    #RMARIN%RASET1%C    #RMARIN%RASET1%CD    #RMARIN%RASET1%CM    #RMARIN%RASET1%I97     #RMARIN%RASET1%J97 !             �@           �                       a                      
      p          p a           p a                                           �@           �                             
                 �@           �                             
                 �@           �                             
                 �@           �                                                �@           �                  !     $                        @                               "                       @                               #            #         @                                   $                    #INDICES %   #NMAX '   #N &            D                                 %                         p          5 � p        r &       5 � p        r &                               
                                  '                     
                                  &           %         @                               (                    
       #RANDUTILS!RANMAR!RASET1 )                                             )     (                    #RANMAR%RASET1%U *   #RANMAR%RASET1%C +   #RANMAR%RASET1%CD ,   #RANMAR%RASET1%CM -   #RANMAR%RASET1%I97 .   #RANMAR%RASET1%J97 /             �@           �                  *     a                      
      p          p a           p a                                           �@           �                  +           
                 �            �                  ,           
                 �            �                  -           
                 �@           �                  .                             �@           �                  /     $             %         @                               0                     
       %         @                                1                     
       %         @                                2                     	                         @                           	     '                    #START_TIME 3   #START 4   #TIME 5   #WRITETIME 6                �                              3                
   1         �   �                       �      4                  #TTIMER_START    1         �   �                      �      5                  #TTIMER_TIME 
   1         �   �                       �      6                  #TTIMER_WRITETIME       �          fn#fn    �   @   J   MPIUTILS !      f       gen@RANDROTATION    f  =       KIND+MPIUTILS !   �  Q       MPISTOP+MPIUTILS %   �  L   a   MPISTOP%MSG+MPIUTILS &   @  R       TTIMER_START+MPIUTILS +   �  T   a   TTIMER_START%THIS+MPIUTILS %   �  Z       TTIMER_TIME+MPIUTILS *   @  T   a   TTIMER_TIME%THIS+MPIUTILS *   �  f       TTIMER_WRITETIME+MPIUTILS /   �  T   a   TTIMER_WRITETIME%THIS+MPIUTILS .   N  L   a   TTIMER_WRITETIME%MSG+MPIUTILS 0   �  @   a   TTIMER_WRITETIME%START+MPIUTILS    �  @       RAND_INST      p       KRAND    �  @       RAND_FEEDBACK    �  V       RANDROTATIONS        $  a   RANDROTATIONS%R     D  @   a   RANDROTATIONS%N    �  V       RANDROTATIOND     �  $  a   RANDROTATIOND%R     �  @   a   RANDROTATIOND%N    >	  W       INITRANDOM    �	  @   a   INITRANDOM%I    �	  @   a   INITRANDOM%I2    
  u       RMARIN (   �
  �   �   RANDUTILS!RMARIN!RASET1     ^  �      RMARIN%RASET1%U       H      RMARIN%RASET1%C !   J  H      RMARIN%RASET1%CD !   �  H      RMARIN%RASET1%CM "   �  H      RMARIN%RASET1%I97 "   "  H      RMARIN%RASET1%J97    j  @   a   RMARIN%IJ    �  @   a   RMARIN%KL    �  f       RANDINDICES $   P  �   a   RANDINDICES%INDICES !     @   a   RANDINDICES%NMAX    D  @   a   RANDINDICES%N    �  m       RANMAR (   �  �   �   RANDUTILS!RANMAR!RASET1     �  �      RANMAR%RASET1%U     i  H      RANMAR%RASET1%C !   �  H      RANMAR%RASET1%CD !   �  H      RANMAR%RASET1%CM "   A  H      RANMAR%RASET1%I97 "   �  H      RANMAR%RASET1%J97    �  P       GAUSSIAN1    !  P       CAUCHY1    q  P       RANDEXP1     �  �       TTIMER+MPIUTILS +   E  H   a   TTIMER%START_TIME+MPIUTILS &   �  Z   a   TTIMER%START+MPIUTILS %   �  Y   a   TTIMER%TIME+MPIUTILS *   @  ^   a   TTIMER%WRITETIME+MPIUTILS 