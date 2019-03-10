podTemplate(
    name: 'jenkins-portcullis',
    label: 'jenkins-portcullis',
    namespace: 'cicd',
    containers: [
      containerTemplate(name: 'cppbuild', image: 'docker.sdlmapleson.net/cppbuild:4', ttyEnabled: true),
      containerTemplate(name: 'gitversion', image: 'docker.sdlmapleson.net/gitversion:5', ttyEnabled: true, command: 'cat'),
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
    stage('Git Version') {
      container('gitversion') {
        sh ("""SEMVER=`semver .` && echo "\$SEMVER" && sed -i "s/AC_INIT(\\[portcullis\\],\\[1.1.2\\]/AC_INIT(\\[portcullis\\],\\[\$SEMVER\\]/" configure.ac""")
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
    stage('Build boost') {
      container('cppbuild') {
        sh "./build_boost.sh"
      }
    }
    stage('Build') {
      container('cppbuild') {
        sh "./autogen.sh && ./configure"
        sh "make -j4 V=1"
      }
    }
    stage('Test') {
      container('cppbuild') {
        echo 'Testing....'
        sh "make -j4 V=1 check"
      }
    }
    stage('Install') {
      container('cppbuild') {
        echo 'Installing....'
        sh "make install"
      }
    }
  }
}
