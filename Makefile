4d:
	rm -rf build && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=4 -DR_NUMERATOR=1 -DR_DENOMINATOR=2 ..

5d:
	rm -rf build && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=5 -DR_NUMERATOR=1 -DR_DENOMINATOR=2 ..
