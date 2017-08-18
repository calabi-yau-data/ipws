dimension=5
r_numerator=1
r_denominator=2
inventory=inventory-tu

suffix=$(dimension)-$(r_numerator)-$(r_denominator)

ansible=ansible -i "$(inventory)"
ansible-playbook=ansible-playbook -i "$(inventory)"

.PHONY: usage
usage:
	@echo "Usage: make ipws pairs deploy copy-pairs find-candidates combine-candidates"

.PHONY: sum
sum:
	$(ansible) all -a "bash -c \"cat /media/data/schoeller/log-5-1-2 | grep unique | awk '{ print $$$$7; }' | paste -sd+ | bc\""

.PHONY: ipws
ipws:
	rm -rf build-$(suffix) && mkdir build-$(suffix) && cd build-$(suffix) && \
	cmake -DCMAKE_BUILD_TYPE=Release -DLINK_STATICALLY=ON \
	      -DDIMENSION=$(dimension) -DR_NUMERATOR=$(r_numerator) \
	      -DR_DENOMINATOR=$(r_denominator) ../.. && \
	make

.PHONY: pairs
pairs:
	mkdir -p data
	build-$(suffix)/ipws-$(suffix) --find-pairs --pairs-out data/pairs-$(suffix) --ws-out data/ws-$(suffix)

.PHONY: deploy
deploy:
	$(ansible-playbook) --extra-vars "suffix=$(suffix)" deploy.yml

.PHONY: copy-pairs
copy-pairs:
	$(ansible-playbook) --extra-vars "suffix=$(suffix)" copy-pairs.yml

.PHONY: find-candidates
find-candidates:
	$(ansible-playbook) --extra-vars "suffix=$(suffix)" find-candidates.yml

.PHONY: combine-candidates
combine-candidates:
	$(ansible-playbook) --extra-vars "suffix=$(suffix)" combine-candidates.yml

.PHONY: ping
ping:
	$(ansible) all -m ping

.PHONY: df
df:
	$(ansible) all -a df

.PHONY: du
du:
	$(ansible) all -a 'du -hs "{{ ws_dir }}"'

.PHONY: count-pairs-used
count-pairs-used:
	$(ansible) all -a "bash -c 'ls /media/data/schoeller/data/ws-$(suffix)-*'" | awk --field-separator=- '{s+=$$6-$$5+1}END{print s}'

.PHONY: count-candidates
count-candidates:
	$(ansible) all -a "bash -c 'cat /media/data/schoeller/log-$(suffix)-*'" | grep unique | awk '{s+=$$7}END{print s}'

# .PHONY: delete-data
# delete-data:
# 	$(ansible) all -a 'rm -r "{{ data_dir }}"'
