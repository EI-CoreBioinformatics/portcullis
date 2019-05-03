podTemplate(
    name: 'jenkins-portcullis',
    label: 'jenkins-portcullis',
    namespace: 'cicd',
    containers: [
      containerTemplate(name: 'cppbuild', image: 'maplesond/cppbuild:latest', ttyEnabled: true),
      containerTemplate(name: 'gitversion', image: 'maplesond/gitversion:latest', ttyEnabled: true, command: 'cat'),
      containerTemplate(name: 'sonarscanner', image: 'docker.sdlmapleson.net/sonarscanner:1', ttyEnabled: true),
      containerTemplate(name: 'docker', image: 'docker', ttyEnabled: true, command: 'cat'),
      containerTemplate(name: 'jnlp', image: 'jenkins/jnlp-slave:3.27-1')
    ],
    volumes: [hostPathVolume(mountPath: '/var/run/docker.sock', hostPath: '/var/run/docker.sock')],
    imagePullSecrets: ['docker.sdlmapleson.net'],
    cloud: 'kubernetes') {
  node('jenkins-portcullis') {
    stage('Git') {
      checkout(scm)
    }
    // Only do semantic versioning on master or develop branch for now
    if(env.BRANCH_NAME == 'master' || env.BRANCH_NAME == 'develop') {
      stage('Git Version') {
        container('gitversion') {
          sh ("""SEMVER=`semver .` && echo "\$SEMVER" && sed -i "s/AC_INIT(\\[portcullis\\],\\[x.y.z\\]/AC_INIT(\\[portcullis\\],\\[\$SEMVER\\]/" configure.ac""")
        }
      }
    }
    // We only have the free sonarqube, so just run on master branch to save time.
    if(env.BRANCH_NAME == 'master') {
      stage('SonarQube') {
        container('sonarscanner') {
          withCredentials([string(credentialsId: 'sonarqube-token', variable: 'TOKEN')]) {
            sh "sonar-scanner -Dsonar.projectKey=portcullis -Dsonar.sources=src,lib,tests,scripts -Dsonar.exclusions=tests/gtest -Dsonar.login=${TOKEN}"
          }
        }
      }
    }    
    stage('Build') {
      container('cppbuild') {
        sh "./autogen.sh"
        sh "./configure"
        sh "make -j4 V=1"
      }
    }
    stage('Test') {
      container('cppbuild') {
        sh "make -j4 V=1 check"
      }
    }
    stage('Install and Package') {
      container('cppbuild') {
        sh "make install"
        sh "make dist"
      }
    }
    if(env.BRANCH_NAME == 'master' || env.BRANCH_NAME == 'develop') {
      stage('Docker') {
        container('docker') {
          def image = docker.build("docker.sdlmapleson.net/portcullis", "--build-arg VERSION=${SEMVER} .")
          image.inside {
            sh "portcullis --help"
            sh "junctools --help"
          }
          docker.withRegistry('https://docker.sdlmapleson.net', 'docker-registry') {
            image.push("${SEMVER}")
            image.push("latest")
            if(env.BRANCH_NAME == 'master') {
              image.push("stable")
            }
          }
        }
      }
    }
  }
}
