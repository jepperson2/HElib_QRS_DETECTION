Vagrant.configure(2) do |config|
  config.vm.box = "ubuntu/trusty64"
  #Using vagrant-disksize plugin to change the size of the disk to 100GB, can adjust to anything >= 40GB
  config.disksize.size = '100GB'
  #config.vm.network "forwarded_port", guest: 5000, host: 5000
  config.vm.provider "virtualbox" do |vb|
    #vb.gui = true
    vb.memory = "7000"
    #The average operation requires 6476MB of RAM for 3 bits!
    #Note that other operations require way less RAM, such as multiplication.
    vb.cpus = 2
    #You can set it to 3 or more to run multiple instances of hqrsd If you have enough RAM.
    vb.name = "hqrsd-virtual-machine"
  end
  config.vm.hostname = "hqrsd-hostname"
  config.vm.provision "shell", inline: <<-SHELL
    sudo apt-get update
    sudo apt-get install -y git g++ m4 perl libboost-all-dev htop 
    sudo apt-get -y autoremove
      
    #Installing GMP
    GMP_V=6.1.2
    wget https://gmplib.org/download/gmp/gmp-$GMP_V.tar.bz2
    tar -xvjf gmp-$GMP_V.tar.bz2
    cd gmp-$GMP_V
    ./configure
    make
    sudo make install
    make check
	cd ..
    rm -fr gmp-$GMP_V*
    unset GMP_V
    
    #Installing NTL
    NTL_V=10.5.0
    wget http://www.shoup.net/ntl/ntl-$NTL_V.tar.gz
    tar -xvzf ntl-$NTL_V.tar.gz
    cd ntl-$NTL_V/src
    ./configure NTL_GMP_LIP=on
    make
    sudo make install
    cd ../..
    rm -fr ntl-$NTL_V*
    unset NTL_V
    
    #Installing HElib
    cd /vagrant
    if [ ! -d "HElib/src" ]; then
        echo "HElib/src not found; Building it again..."
        #rm -fr HElib
        git clone https://github.com/shaih/HElib.git
        cd HElib/src
        make
        #make check
		#1 of the check of HElib fails for some reason
        #make test
        cd ../..
    fi
   
    #Building hqrsd
    make hqrsd
    echo "hqrsd was successfully built ! Run it with: cd /vagrant && ./hqrsd"
  SHELL
end

