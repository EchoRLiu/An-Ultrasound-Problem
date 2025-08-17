# An Ultrasound Problem

**Abstract**

Most of the time, data set comes with noise, covering the true information we hope to extract. For instance, when analyzing ultrasound information, the real signal is always buried somewhere. Here, we present two steps to denoise the data and detect the locations of certain moving object. The first step is to take the average of Fourier transform of the data through realizations. Doing this in frequency domain instead of spatial domain is because the movement of the object in spatial domain would obscure the result, but instead, in frequency domain, since the frequency of the object is fixed, through averaging, we would not add noise to the data. This would also give us the center frequency of the object, which would be used in the second step to clean the data further. The second step is to use Gaussian filter to filter the data in the frequency domain, and get the location of the object through the inverse Fourier transform of the filtered result.

**Introduction and Overview**

Here, we use the following example problem to illustrate how Fourier transform can be used to denoise the data and locate the object of interest.

**Problem Description**

In this example problem, our dog swallowed a marble by accident. Through ultrasound around the area we suspect the marble would be, we obtain 20 realizations of the spatial domain data of the marble locations. But as always, other signals in the environment and the movement of our dog would generate highly noise in our data. In addition, through realization, the marble would also be moving, which adds difficulty in denoising our data.

**General Approach**

Considering the movement of the marble through realizations, we can not simply denoise the data by averaging the spatial domain data (we assume the noise would be white noise, i.e. they are independent normal distribution, and if we average them, the mean should be zero). Thus, we first obtain the Fourier transform of the data, which represents all frequencies we have in the data. We then average the frequency domain data, to extract the center frequency, with which we can further denoise the data. Using 3-D Gaussian filter around the center frequency, we are able to clean up the data. The last step would be to use inverse Fourier transform to obtain the location of the marble.

**Theoretical Background**

As we learned from our [textbook](https://faculty.washington.edu/kutz/582.pdf), Fourier introduced the concept of representing a given function f(x) by a trigonometric series of sines and cosines:

    $f(x) = \frac{a_0}{2} + \sum_{i=1}^\infty \left(a_n\cos{nx} + b_n\sin{nx}\right) \quad x \in (-\pi,\pi]$

And further, $F(k)$, the Fourier transform of a function $f(x)$:

    F(k) = \frac{1}{\sqrt{2\pi}} \int_{-\inf}^{\inf} e^{-ikx}f(x)dx

is able to represent the frequency information in $f(x)$. Consider $f(x)$ as a combination of some signals. If we look at the Fourier transform of every individual signal in $f(x)$, it would only have a spike value around the frequency of this specific signal, meaning $F(k)$ would have a relatively large value around $k$ if k is the frequency of the signal. Now as for the Fourier transform of $f(x)$, which represents the combination of all the signals, it would gives us the information of all the frequencies we have in $f(x)$.

An important aspect of denoising mentioned in the [textbook](https://faculty.washington.edu/kutz/582.pdf), is frequency filtering. One simple type of filter is Gaussian filter:

    F(k) = exp(-\tau (k-k_0)^2)

which would basically eliminate frequency far away from $k_0$. If $k_0$ is the frequency of the object we are interested in, the inverse Fourier transform of the filtered result would give us very nicely denoised data. Otherwise, the inverse Fourier transform of the filtered result would be near $0$, since the filtered result is basically white noise. The main idea is that by using Gaussian filter around the frequency of interest, we are able to clear the noise in the data, and obtain the important information we hope to extract.

**Algorithm Implementation and Development**

The general approach is implemented in the following way.

The ultrasound data is loaded into $Undata$ and the spatial domain length $L$ and Fourier modes $n$ are defined for discretization. The spatial and spectral coordinates are also created. Since $fft$ command in Matlab returns a shifted version, to see the correct mathematical result, we need to use $fftshift$. Here, to save future multiple computation, we also store $ks$, as the correct version of $k$, meaning it's shifted back to the correct version.
    
Through visualisation of 20 realizations of data, we can clearly see the data is very noisy.  
    
The data is then transformed into frequency domain using Algorithm 1.

    FOR j = 1:20
        Extract realization $j$ from \texttt{Undata} and reshape into (n,n,n) as $Un_j$
        Get the Fourier transform version of $Un_j$ as $Unt_j$
        shift $Unt_j$ to the right mathematical version for visualisation
        Visualize the data
    END
    
To clean the white noise a bit, we average the Fourier transformed result for each realization using Algorithm 2. We then visualize the averaged result, and try to extract the center frequency.
    
    Initialize $Untave$ to store the averaged result of Fourier transforms
    FOR j = 1:20
        Extract realization $j$ from \texttt{Undata} and reshape into (n,n,n) as $Un_j$
        Get the Fourier transform version of $Un_j$ as $Unt_j$
        Add $Unt_j$ to $Untave$
    END
    shift $Untave$ to the right mathematical version for visualisation and divide by $20$ to get the correct average as $Untave$.
    Visualize the data
    Search for the coordinates of spike value in the averaged Fourier transforms to find center frequency.
    
We then try to use 3D Gaussian filter to further clean the data using Algotithm 3. We also can obtain the locations of the marble along the $20$ realizations.
    
    Define the 3D Gaussian filter around $\omega_x = 1.89, \omega_y = -1.05, \omega_z = 0.0$.
    Initialize $marble_pos$ to store the spatial coordinates of the marble.
    FOR j = 1:20$
        Extract realization $j$ from \texttt{Undata} and reshape into (n,n,n) as $Un_j$
        Get the Fourier transform version of $Un_j$ as $Unt_j$
        Filter the data: $Unft_j$ = filter*$Unt_j$
        Obtain the inverse Fourier Transform of $Unft_j$ as $Unf_j$
        Visualize the cleaned data
        Obtain the coordinates of the spike value in $Unf_j$ as the location of the marble in realization j.
    END

**Computational Results**

At first, the data is very noisy. Figure 1 is an example of the realizations in the spatial domain, and Figure 2 is one in the frequency domain.
    
![figure 1](https://github.com/EchoRLiu/An-Ultrasound-Problem/blob/master/noisydata1.jpg)
    
![figure 2](https://github.com/EchoRLiu/An-Ultrasound-Problem/blob/master/noisydata2.jpg)
    
After we average the frequency domain data, we can see in Figure 3, the result seems less noisy compared to Figure 2. We can also notice the center frequency is near $\omega_x = 1.89, \omega_y = -1.05, \omega_z = 0.0$.
    
![figure 3](https://github.com/EchoRLiu/An-Ultrasound-Problem/blob/master/averagedWN2.jpg)

To compare the results between averaged spatial and frequency domains, see Figure 4.

![figure 4](https://github.com/EchoRLiu/An-Ultrasound-Problem/blob/master/averagedWN1.jpg)
    
After filtering, we can clearly see the marble's moving path as in Figure 5. A simple path is also presented here in Figure 6.
    
![figure 5](https://github.com/EchoRLiu/An-Ultrasound-Problem/blob/master/marbles.jpg)
![figure 6](https://github.com/EchoRLiu/An-Ultrasound-Problem/blob/master/marblepath.jpg)
    
Finally, in the 20th realization, the position of marble we obtained is $x = -5.6250, y = 4.2188, z = -6.0938$, which is the location where the intense acoustic wave should be focused.

**Summary and Conclusions**

Through this example problem, we can see how useful Fourier transform is for signal analysis and noise attenuation. Fourier transform is an amazing tool that makes it possible for us to separate mixed signals and analyze individually.
