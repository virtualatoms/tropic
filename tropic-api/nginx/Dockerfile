FROM nginx:1.21-alpine

# Remove default configuration
RUN rm /etc/nginx/conf.d/default.conf

# Copy our custom configuration
COPY default.conf /etc/nginx/conf.d/default.conf

# Copy SSL certificates (if using HTTPS)
# COPY ssl/ /etc/nginx/ssl/

EXPOSE 80
# EXPOSE 443 # Uncomment when using HTTPS
