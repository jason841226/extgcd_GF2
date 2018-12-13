import random
f=open('random_input.txt','w');
for i in range(0,32):
	f.write(str(random.randint(0,1)));
	if i%8==7:
		f.write(' ');
