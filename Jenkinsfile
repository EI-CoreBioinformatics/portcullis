podTemplate(
    name: 'jenkins-portcullis',
    label: 'jenkins-portcullis',
    namespace: 'cicd',
    containers: [
      containerTemplate(name: 'cppbuild', image: 'maplesond/cppbuild:latest', ttyEnabled: true),
      containerTemplate(name: 'gitversion', image: 'maplesond/gitversion:latest', ttyEnabled: true, command: 'cat'),
      containerTemplate(name: 'sonarscanner', image: 'docker.sdlmapleson.net/sonarscanner:1', ttyEnabled: true),
      containerTemplate(name: 'githubrelease', image: 'maplesond/githubrelease:latest', ttyEnabled: true, command: 'cat'),
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
          sh ("""SEMVER=`semver .` && echo "\$SEMVER" > version && sed -i "s/AC_INIT(\\[portcullis\\],\\[x.y.z\\]/AC_INIT(\\[portcullis\\],\\[\$SEMVER\\]/" configure.ac""")
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
        if(env.BRANCH_NAME == 'master') {
          sh "echo 'TODO send release to github'"
        }        
      }
    }
    if(env.BRANCH_NAME == 'master' || env.BRANCH_NAME == 'develop') {
      stage('Docker') {
        container('docker') {
          // Get the version we previously saved to disk
          def versionFile = readFile "version"
          SEMVER = versionFile.split('\n')[0]
          // Build the image
          def image = docker.build("docker.sdlmapleson.net/portcullis:${SEMVER}", "--build-arg VERSION=${SEMVER} .")
          image.inside {
            sh "portcullis --help || true"
            sh "junctools --help || true"
          }
          // Push to local registry
          docker.withRegistry('https://docker.sdlmapleson.net', 'docker-registry') {
            image.push()
            image.push("latest")
            if(env.BRANCH_NAME == 'master') {
              withCredentials([usernamePassword(credentialsId: 'dockerhub', usernameVariable: 'DOCKERHUB_USER', passwordVariable: 'DOCKERHUB_PASS')]) {
                sh "docker push docker.sdlmapleson.net/portcullis:stable"
                sh "docker tag docker.sdlmapleson.net/portcullis:${SEMVER} maplesond/portcullis:${SEMVER}"
                sh "docker tag docker.sdlmapleson.net/portcullis:${SEMVER} maplesond/portcullis:latest"
                sh "docker tag docker.sdlmapleson.net/portcullis:${SEMVER} maplesond/portcullis:stable"
                sh "docker login -u ${DOCKERHUB_USER} -p ${DOCKERHUB_PASS}"
                sh "docker push maplesond/portcullis:${SEMVER}"
                sh "docker push maplesond/portcullis:latest"
                sh "docker push maplesond/portcullis:stable"
              }
            }
          }
        }
      }
    }
    if(env.BRANCH_NAME == 'master') {
      stage('Github Release') {
        withCredentials([string(credentialsId: 'github-maplesond-portcullis', variable: 'GITHUB_TOKEN')]) {
          sshAgent('gitea') {
            def versionFile = readFile "version"
            SEMVER = versionFile.split('\n')[0]
            sh "git tag Release-${SEMVER} && git push --tags"
            sh "git remote add github-maplesond https://github.com/maplesond/portcullis.git && git push github-maplesond master && git push github-maplesond master --tags"
            // Ignore these for now... not sure if they are being used.
            //sh "git remote github-ei https://github.com/EI-CoreBioinformatics/portcullis.git"
            //sh "git remote github-tgac https://github.com/TGAC/portcullis.git"
            GITHUB_USER=maplesond
            GITHUB_REPO=portcullis
            sh "DESCRIPTION=`git log -1 | tail -n +4`"
            sh "github-release releases --tag ${SEMVER}"
            sh "github-release upload --tag ${SEMVER} --name ${GITHUB_REPO}-${SEMVER}.tar.gz --file ${GITHUB_REPO}-${SEMVER}.tar.gz --label 'source code distributable' --description ${DESCRIPTION}"
          }
        }
      }
      // This is probably going to be hard but it should be possible to at least start the process in an automatic way.
      //stage('Brew update') {
      //}
      //stage('Conda update') {
      //}
    }
  }
}
