function y = convolution(element1, element2)

foo = fft(element1);
foo = foo .* fft(element2);
y = ifft(foo);