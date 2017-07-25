dimension=6
r_numerator=3
r_denominator=2

suffix=$(dimension)-$(r_numerator)-$(r_denominator)

ansible=ansible -i inventory
ansible-playbook=ansible-playbook -i inventory

.PHONY: usage
usage:
	@echo Usage: make ipws pairs deploy find-ws

.PHONY: ipws
ipws:
	rm -rf build-$(suffix) && mkdir build-$(suffix) && cd build-$(suffix) && \
	cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON \
	      -DDIMENSION=$(dimension) -DR_NUMERATOR=$(r_numerator) \
	      -DR_DENOMINATOR=$(r_denominator) ../.. && \
	make

pairs:
	mkdir -p data
	build-$(suffix)/ipws-$(suffix) --find-pairs --pairs-out data/pairs-$(suffix) --ws-out data/ws-$(suffix)

deploy:
	$(ansible-playbook) --extra-vars "suffix=$(suffix)" deploy.yml

find-ws:
	$(ansible-playbook) --extra-vars "suffix=$(suffix)" find-ws.yml

ping:
	$(ansible) all -m ping

df:
	$(ansible) all -a df

du:
	$(ansible) all -a 'du "{{ work_dir }}"'
