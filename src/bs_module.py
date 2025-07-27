import jax.numpy as jnp
from jax import grad, vmap
from jax.scipy.stats import norm as jnorm

class BSModel:

    def __init__(self, r: float, q: float = 0.0):
        self.r = r
        self.q = q
        # Vektorisering av pris- og greeks-funksjoner
        self.price_vec = vmap(self.price, in_axes=(0,0,0,0,None))
        self.delta_vec = vmap(self.delta, in_axes=(0,0,0,0,None))
        self.vega_vec = vmap(self.vega, in_axes=(0,0,0,0,None))

    def price(self, S, K, T, sigma, option: str = 'call'):
        d1 = (jnp.log(S/K) + (self.r - self.q + 0.5 * sigma**2) * T) / (sigma * jnp.sqrt(T))
        d2 = d1 - sigma * jnp.sqrt(T)
        dfS = jnp.exp(-self.q * T)
        dfK = jnp.exp(-self.r * T)
        if option == 'call':
            return S * dfS * jnorm.cdf(d1) - K * dfK * jnorm.cdf(d2)
        elif option == 'put':
            return K * dfK * jnorm.cdf(-d2) - S * dfS * jnorm.cdf(-d1)
        else:
            raise ValueError("option must be 'call' or 'put'")

    def delta(self, S, K, T, sigma, option: str = 'call'):
        d1 = (jnp.log(S/K) + (self.r - self.q + 0.5 * sigma**2) * T) / (sigma * jnp.sqrt(T))
        dfS = jnp.exp(-self.q * T)
        return dfS * (jnorm.cdf(d1) if option=='call' else jnorm.cdf(d1)-1)

    def vega(self, S, K, T, sigma, option: str = 'call'):
        d1 = (jnp.log(S/K) + (self.r - self.q + 0.5 * sigma**2) * T) / (sigma * jnp.sqrt(T))
        dfS = jnp.exp(-self.q * T)
        return S * dfS * jnorm.pdf(d1) * jnp.sqrt(T)

    def implied_vol(self, market_price, S, K, T, option='call',
                     sigma0=0.2, tol=1e-6, maxiter=50):
        """
        Newton-Raphson for å finne implied vol.
        Støtter skalar eller vektor via explicitt loop.
        """
        # Hvis arrays: iterér
        if hasattr(market_price, 'shape') and market_price.shape != ():
            ivs = []
            for mp, s, k, t in zip(market_price, S, K, T):
                ivs.append(self.implied_vol(mp, s, k, t, option, sigma0, tol, maxiter))
            return jnp.array(ivs)
        # Skalar-case
        def price_fn(sig):
            return self.price(S, K, T, sig, option)
        grad_fn = grad(price_fn)
        sigma = sigma0
        for i in range(maxiter):
            f = price_fn(sigma) - market_price
            v = grad_fn(sigma)
            update = f / v
            sigma -= update
            if jnp.abs(update) < tol:
                break
        return jnp.abs(sigma)























