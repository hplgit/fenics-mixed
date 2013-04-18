      subroutine set_ic_and_dt(n, s0, dt)
Cf2py intent(c) set_ic_and_dt
      integer n
      real*8 s0(0:n-1), dt
Cf2py intent(c) n, s0, dt
      return
      end

      subroutine set_const_ic_and_dt(n, s0, dt)
Cf2py intent(c) set_const_ic_and_dt
      integer n
      real*8 s0, dt
Cf2py intent(c) n, s0, dt
      return
      end

      subroutine set_u(n, u)
Cf2py intent(c) set_u
      integer n
      real*8 u(0:n-1)
Cf2py intent(c) n, u
      return
      end

      subroutine advance_one_timestep()
Cf2py intent(c) advance_one_timestep
      return
      end

      subroutine get_s(n, s)
Cf2py intent(c) get_s
      integer n
Cf2py intent(c) n
      real*8 s(0:n-1)
Cf2py intent(c, in, out) s
      return
      end
