3-1-2:
	rm -rf build-3-1-2
	mkdir build-3-1-2
	cd build-3-1-2 && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=3 -DR_NUMERATOR=1 -DR_DENOMINATOR=2 ..

4-1-1:
	rm -rf build-4-1-1
	mkdir build-4-1-1
	cd build-4-1-1 && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=4 -DR_NUMERATOR=1 -DR_DENOMINATOR=1 ..

4-1-2:
	rm -rf build-4-1-2
	mkdir build-4-1-2
	cd build-4-1-2 && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=4 -DR_NUMERATOR=1 -DR_DENOMINATOR=2 ..

5-1-1:
	rm -rf build-5-1-1
	mkdir build-5-1-1
	cd build-5-1-1 && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=5 -DR_NUMERATOR=1 -DR_DENOMINATOR=1 ..

5-1-2:
	rm -rf build-5-1-2
	mkdir build-5-1-2
	cd build-5-1-2 && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=5 -DR_NUMERATOR=1 -DR_DENOMINATOR=2 ..

6-3-2:
	rm -rf build-6-3-2
	mkdir build-6-3-2
	cd build-6-3-2 && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=6 -DR_NUMERATOR=3 -DR_DENOMINATOR=2 ..

6-1-1:
	rm -rf build-6-1-1
	mkdir build-6-1-1
	cd build-6-1-1 && cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON -DDIMENSION=6 -DR_NUMERATOR=1 -DR_DENOMINATOR=1 ..

clean:
	rm -rf build-*
