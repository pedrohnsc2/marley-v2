#!/usr/bin/env bash
# =============================================================================
# Marley VPS Setup Script
#
# Prepares a fresh Ubuntu 22.04/24.04 VPS for running Marley.
# Tested on: DigitalOcean, Hetzner, Oracle Cloud Free Tier.
#
# Usage:
#   ssh root@your-server
#   curl -fsSL https://raw.githubusercontent.com/YOUR_REPO/main/infra/setup-vps.sh | bash
#   # or: scp this file to server and run it
# =============================================================================

set -euo pipefail

echo "=== Marley VPS Setup ==="
echo "Date: $(date)"
echo "Hostname: $(hostname)"

# ---------------------------------------------------------------------------
# 1. System packages
# ---------------------------------------------------------------------------
echo "--- Installing system packages ---"
apt-get update
apt-get install -y \
    curl \
    git \
    ufw \
    fail2ban \
    htop \
    unattended-upgrades

# ---------------------------------------------------------------------------
# 2. Docker
# ---------------------------------------------------------------------------
if ! command -v docker &> /dev/null; then
    echo "--- Installing Docker ---"
    curl -fsSL https://get.docker.com | bash
    systemctl enable docker
    systemctl start docker
fi

# Docker Compose V2 plugin (comes with modern Docker install)
docker compose version || {
    echo "ERROR: docker compose not available"
    exit 1
}

# ---------------------------------------------------------------------------
# 3. Non-root deploy user
# ---------------------------------------------------------------------------
if ! id marley &> /dev/null; then
    echo "--- Creating marley user ---"
    useradd -m -s /bin/bash -G docker marley
    echo "NOTE: Set up SSH key for marley user:"
    echo "  mkdir -p /home/marley/.ssh"
    echo "  cp ~/.ssh/authorized_keys /home/marley/.ssh/"
    echo "  chown -R marley:marley /home/marley/.ssh"
fi

# ---------------------------------------------------------------------------
# 4. Firewall
# ---------------------------------------------------------------------------
echo "--- Configuring firewall ---"
ufw --force reset
ufw default deny incoming
ufw default allow outgoing
ufw allow ssh
ufw allow 80/tcp
ufw allow 443/tcp
ufw --force enable
ufw status

# ---------------------------------------------------------------------------
# 5. Fail2ban
# ---------------------------------------------------------------------------
echo "--- Configuring fail2ban ---"
systemctl enable fail2ban
systemctl start fail2ban

# ---------------------------------------------------------------------------
# 6. Swap (for small VPS — 2GB swap)
# ---------------------------------------------------------------------------
if [ ! -f /swapfile ]; then
    echo "--- Creating 2GB swap ---"
    fallocate -l 2G /swapfile
    chmod 600 /swapfile
    mkswap /swapfile
    swapon /swapfile
    echo '/swapfile swap swap defaults 0 0' >> /etc/fstab
fi

# ---------------------------------------------------------------------------
# 7. Application directory
# ---------------------------------------------------------------------------
echo "--- Setting up application directory ---"
mkdir -p /opt/marley/infra
chown -R marley:marley /opt/marley

# ---------------------------------------------------------------------------
# 8. Auto-updates
# ---------------------------------------------------------------------------
echo "--- Enabling unattended security upgrades ---"
dpkg-reconfigure -f noninteractive unattended-upgrades

# ---------------------------------------------------------------------------
# 9. Docker log rotation
# ---------------------------------------------------------------------------
echo "--- Configuring Docker log rotation ---"
cat > /etc/docker/daemon.json << 'EOF'
{
  "log-driver": "json-file",
  "log-opts": {
    "max-size": "10m",
    "max-file": "3"
  }
}
EOF
systemctl restart docker

echo ""
echo "=== VPS Setup Complete ==="
echo ""
echo "Next steps:"
echo "  1. Copy .env.production to /opt/marley/"
echo "  2. Set up SSH key for GitHub Actions deploy"
echo "  3. Push to main to trigger first deploy"
echo "  4. Or manually: cd /opt/marley && docker compose up -d"
echo ""
echo "Recommended VPS (South America optimized):"
echo "  - Hetzner CPX21: 3 vCPU, 4GB RAM, ~EUR 8/mo (Ashburn DC, low latency to BR)"
echo "  - Oracle Cloud Free Tier: ARM 4 vCPU, 24GB RAM, 200GB (Sao Paulo region!)"
echo "  - DigitalOcean: 2 vCPU, 4GB RAM, ~USD 24/mo (no BR region)"
