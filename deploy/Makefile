dimension=6
r_numerator=3
r_denominator=2

suffix=$(dimension)-$(r_numerator)-$(r_denominator)

ansible=ansible -i inventory
ansible-playbook=ansible-playbook -i inventory

.PHONY: ipws
ipws:
	rm -rf build && mkdir build && cd build && \
	cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON \
	      -DDIMENSION=$(dimension) -DR_NUMERATOR=$(r_numerator) \
	      -DR_DENOMINATOR=$(r_denominator) ../.. && \
	make

pairs:
	mkdir -p data
	build/ipws --find-pairs --pairs-out data/pairs-$(suffix) --ws-out data/ws-$(suffix)

.PHONY: count-pairs
count-pairs:
	$(eval pair_count=`build/ipws --print-count --pairs-in data/$(pairs)`)

deploy: count-pairs
	$(ansible-playbook) --extra-vars "suffix=$(suffix)" deploy.yml

ping:
	$(ansible) all -m ping

df:
	$(ansible) all -a df
