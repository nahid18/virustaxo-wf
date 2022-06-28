FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:02ab-main
RUN python3 -m pip install biopython numpy pandas python-dateutil pytz six tqdm
COPY wf /root/wf

ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN  sed -i 's/latch/wf/g' flytekit.config
RUN python3 -m pip install --upgrade latch
WORKDIR /root
