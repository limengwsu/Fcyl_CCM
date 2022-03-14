 function F = visco(Sal, Pp, TC)
 %calculate viscosity of seawater - linked to diffusion coeficient code
 % Valid for 0<T<30 and 0<S<36, Calculated viscosity is in centipoise.
 
 F =  1.7910 - TC*(6.144D-02 - TC*(1.4510D-03 - TC*1.6826D-05))...
       - 1.5290D-04*Pp + 8.3885D-08*Pp^2 + 2.4727D-03*Sal...
       + (6.0574D-06*Pp - 2.6760D-09*Pp^2)*TC + (TC*(4.8429D-05...
       - TC*(4.7172D-06 - TC*7.5986D-08)))*Sal;
   
   return