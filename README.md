# **Analysis of photometric data of an observed star**

## **ENG**
This code was used for the photometric reduction of an observed star and calculation of its light curve.
The sequence of operations performed is as follows (The observed star used as an example is NY Vir star.):

1. Data reduction was performed for NY Vir images (Master Bias and Master Flat).
2. Additional image cleaning and trimming were performed for NY Vir images.
3. Trimmed images were aligned.
4. Photometry was performed by selecting the pixel position of NY Vir and a nearby star that is similar in brightness and position.
5. The light curve of NY Vir was drawn by measuring its differential magnitude in the Fark Işık band, and an appropriate fit model was drawn to determine the minimum time.
6. The results of the light curve image and the drawn fit model were recorded.

## **TR**
Bu kod, NY Vir star'ın fotometrik indirgemesi yapılıp ışık eğrisi hesaplanması için kullanılmıştır. İşlemler şöyle yapılmıştır:

1. NY Vir görüntüleri için veri indirgemesi yapıldı (master bias ve master flat).
2.NY Vir görüntüleri için ekstra görüntü temizleme ve trim yapıldı.
3. Trimlenmiş görüntüler hizalandı.
4. NY Vir ve ona hem parlaklık hem de konum olarak yakın bir yıldızın piksel konumu seçilerek fotometri yapıldı.
5. Fotometri yapılan NY Vir'in fark ışık ölçümü yapılarak ışık eğrisi çizildi ve minimum zamanı belirlenmesi için uygun bir fit modeli çizildi.
6. Işık eğrisi görüntüsü ve çekilen fit modelinin sonuçları kaydedildi.

Bu kod, veri indirgeme, görüntü işleme ve fotometri analizi yaparak, astronomi araştırmalarında kullanılan standart bir yöntem olan "light curve" analizini gerçekleştirir.
Işık eğrisi, gözlemlenen bir yıldızın ışık şiddetinin zamana göre değişimini gösteren bir grafiktir.
Bu analiz, gözlemlenen yıldızın parlaklık değişimlerinin incelenmesi ve yıldızın fiziksel özellikleri hakkında bilgi edinilmesi için önemlidir.

### **RESULT / SONUÇ**:
    
    ```
    [[Model]]
    (Model(gaussian, prefix='g1_') + Model(polynomial, prefix='p_'))
    [[Fit Statistics]]
    # fitting method = leastsq
    # function evals = 635
    # data points = 77
    # variables = 5
    chi-square = 0.01715809
    reduced chi-square = 2.3831e-04
    Akaike info crit = -637.500023
    Bayesian info crit = -625.780996
    [[Variables]]
    p_c0: 1.15515834 +/- 0.14938206 (12.93%) (init = 1.05553)
    p_c1: 0.66941431 +/- 0.29920806 (44.70%) (init = 0.4445682)
    g1_amplitude: -0.00457827 +/- 4.4681e-05 (0.98%) (init = 1)
    g1_center: 0.49927357 +/- 1.8251e-05 (0.00%) (init = 0.5)
    g1_sigma: 0.00235736 +/- 2.0801e-05 (0.88%) (init = 0.01)
    g1_fwhm: 0.00555116 +/- 4.8983e-05 (0.88%) == '2.3548200*g1_sigma'
    g1_height: -0.8 (fixed)
    [[Correlations]] (unreported correlations are < 0.500)
    C(p_c0, p_c1) = -1.000
    C(g1_amplitude, g1_sigma) = -0.737
    ```

![fotometric_fit](/NYVir.png)