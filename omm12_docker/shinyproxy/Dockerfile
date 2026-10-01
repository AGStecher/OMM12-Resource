# ShinyProxy: routes each visitor to their own isolated container running
# the omm12 image, so multiple users can use the app in parallel without
# sharing one R session. Same pattern as daisybio/namco's shinyproxy image.
# (namco's openjdk:11-jre base image has been removed from Docker Hub;
#  eclipse-temurin is its official successor.)
FROM eclipse-temurin:11-jre

RUN mkdir -p /opt/shinyproxy/
ADD https://www.shinyproxy.io/downloads/shinyproxy-2.4.0.jar /opt/shinyproxy/shinyproxy.jar
# application.yml (local config, git-ignored) is used if present; otherwise the
# template from the repository. With authentication: none the admin login in it
# is not used, so the template works as is.
COPY application*.yml /opt/shinyproxy/
RUN cd /opt/shinyproxy && \
    if [ ! -f application.yml ]; then mv application.example.yml application.yml; \
    else rm -f application.example.yml; fi

WORKDIR /opt/shinyproxy/
CMD ["java", "-jar", "/opt/shinyproxy/shinyproxy.jar"]
