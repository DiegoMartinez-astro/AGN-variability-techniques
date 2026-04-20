import numpy as np

# CLASS LIGHT CURVE
class LC:
    def __init__(self, t, y, yerr, label="AGN Data",color='mediumseagreen'):
        """Basic light curve parameters"""
        self.t = np.array(t)
        self.y = np.array(y)
        self.yerr = np.array(yerr)
        self.label = label
        self.color = color

    """WEIGHTED BINNING FUNCTION"""

    def binning(self, bin_size=1):

       # Define limits
       t_min = np.min(self.t)
       t_max = np.max(self.t)

       # Define new time array using bin size
       t_new = np.arange(t_min,t_max+bin_size,bin_size)

       t_bin = []
       y_bin = []
       yerr_bin = []

       for i in range(len(t_new)-1):
           low_lim = t_new[i]
           upp_lim = t_new[i+1]

           mask = (self.t <  upp_lim) & (self.t >= low_lim) #bin mask
           
            #collect elements within current bin
           t_mask = self.t[mask]
           y_mask = self.y[mask]
           yerr_mask = self.yerr[mask]

           if np.any(mask):

               # Weight operations using the magnitude/flux errors
               w = 1/(yerr_mask**2)

               # Define mean value for the bin
               y_mean = (np.sum(y_mask*w))/np.sum(w)

               # error
               yerr_mean = np.sqrt(1/(np.sum(w)))

               # new time
               t_mean = np.mean(t_mask)

               # add bin values 
               t_bin.append(t_mean)
               y_bin.append(y_mean)
               yerr_bin.append(yerr_mean)

        # transform everything to numpy array

       t_bin = np.array(t_bin)
       y_bin = np.array(y_bin)
       yerr_bin = np.array(yerr_bin)

        # save result as a new light curve object
       result = LC(t_bin, y_bin, yerr_bin, color = self.color, label=self.label+''+str(bin_size)+' day binned')

       return result

    """UNDER SAMPLING FUNCTION"""

    def decadence(self, target_cadence):
        t_dec, y_dec, yerr_dec = [], [], []
    
        # define time bins
        bins = np.arange(np.min(self.t), np.max(self.t) + target_cadence, target_cadence)

        for i in range(len(bins)-1):
            # index all fluxes whithin a bin
            idx = np.where((self.t >= bins[i]) & (self.t < bins[i+1]))[0]
        
            # Avoid empty bins
            
            if len(idx) > 0:
                # choose a random element from the bin
                chosen_index = np.random.choice(idx)
            
            # Save random element
                t_dec.append(self.t[chosen_index])
                y_dec.append(self.y[chosen_index])
                yerr_dec.append(self.yerr[chosen_index])
        
        result = LC(t_dec, y_dec, yerr_dec, color = self.color, label=self.label+''+str(target_cadence)+' days cadence')
        return result
    
    
    """SIGMA CLIPPING"""
        
    def sigma_clip(self, sigma=3):

        mean = np.mean(self.y)
        std = np.std(self.y)

        # Limits
        lower_bound = mean - sigma * std
        upper_bound = mean + sigma * std

        # filter
        mask = (self.y >= lower_bound) & (self.y <= upper_bound)

        # removed points
        
        n_removed = len(self.y) - np.sum(mask)
        
        print(f"Sigma-clipping: {n_removed} points of {len(self.y)} have been removed.")

        result = LC(self.t[mask], self.y[mask], self.yerr[mask], label=self.label, color=self.color)
        

        return result  #return cleaned lightcurve

    
    """ J-INDEX AS DEFINED IN Ma+2024 """

    def j_index(self):
        
        n=len(self.y)
        
        mean = np.mean(self.y)

        #delta value
        delta = ((self.y - mean) / self.yerr) * np.sqrt(n / (n - 1))

        #weights
        weight = (1 + (delta / 2)**2)**(-1)

        #numerator of j
        term = weight * np.sign(delta**2 - 1) * np.sqrt(np.abs(delta**2 - 1))
        
        #final result
        j = np.sum(term) / np.sum(weight)

        return j


    """ SMOOTHNESS PARAMETER AS DEFINED IN D. Martinez Collipal & S. Panda 2026 """
    
    def s(self):

        t = self.t

        # get mean flux
        mean = np.mean(self.y)

        value=0

        N = len(self.t)

        n=0

        #start s calculation
        
        for i in range(0,len(self.y)-2):
            
            
            # dt interval in days to apply MA+2024 condition
            
            dt= self.t[i+2]-self.t[i]

            # Calculate smoothness in regions without gaps
            if dt < 20:

                #numerators 
                
                df_1=self.y[i+1]-self.y[i]
                df_2=self.y[i+2]-self.y[i+1]

                #denominators 
                
                dt_1=t[i+1]-t[i]
                dt_2=t[i+2]-t[i+1]

                #delta f_i

                delta_1 = (1/mean)*(df_1/dt_1)

                delta_2 = (1/mean)*(df_2/dt_2)

                #separate sum in two components
                
                A = np.abs((delta_2-delta_1)/(t[i+2]-t[i]))

                B = (((df_2)**2 + (df_1)**2 )**0.5)/mean

                #final element

                value+= A*B
                n+=1
            else:
                continue

        value_n = value/N
        return value_n
        

    def __repr__(self):
        return f"LightCurve(points={len(self.t)}, label='{self.label}')"
